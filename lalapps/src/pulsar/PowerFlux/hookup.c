#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <regex.h>

#include <lal/LALInitBarycenter.h>
#include <lal/DetResponse.h>
#include <lal/Velocity.h>
#include <lal/DetectorSite.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "global.h"
#include "cmdline.h"
#include "grid.h"
#include "intervals.h"
#include "polarization.h"

#ifndef PATH_MAX
/* just in case it is not defined */
#define PATH_MAX PATH_MAX
#endif

extern FILE *LOG, *FILE_LOG;
extern char *earth_ephemeris;
extern char *sun_ephemeris;
extern struct gengetopt_args_info args_info;
extern char *output_dir;

INTERVAL_SET *segment_list=NULL;
INTERVAL_SET *veto_segment_list=NULL;

regex_t write_dat, write_png;

EphemerisData ephemeris;
LALDetector detector;

void init_hookup(void)
{
if(regcomp(&write_dat, args_info.write_dat_arg, REG_EXTENDED | REG_NOSUB)){
	fprintf(stderr,"Cannot compile \"--write-dat=%s\"\n", args_info.write_dat_arg);
	exit(-1);
	}
if(regcomp(&write_png, args_info.write_png_arg, REG_EXTENDED | REG_NOSUB)){
	fprintf(stderr,"Cannot compile \"--write-dat=%s\"\n", args_info.write_dat_arg);
	exit(-1);
	}
if(args_info.segments_file_given){
	segment_list=new_interval_set();
	add_intervals_from_file(segment_list, args_info.segments_file_arg);
	fprintf(LOG, "Read %d intervals from file: %s\n", segment_list->free, args_info.segments_file_arg);
	}
if(args_info.veto_segments_file_given){
	veto_segment_list=new_interval_set();
	add_intervals_from_file(veto_segment_list, args_info.veto_segments_file_arg);
	fprintf(LOG, "Read %d veto intervals from file: %s\n", veto_segment_list->free, args_info.veto_segments_file_arg);
	}
}

int get_power_range(char *filename, long startbin, long count, float *data, INT64 *gps)
{
FILE *fin;
char s[PATH_MAX];
char *p;
errno=0;
fin=fopen(filename,"r");
if(errno==ENOENT)return -1; /* no such file */
if(fin==NULL){
	fprintf(stderr,"Error reading band %ld-%ld bins from file \"%s\":",startbin,startbin+count,filename);
	perror("");
	return -1;	
	}
memset(s,0,PATH_MAX);
while(strncmp(s,"binary",6)){
	fgets(s,19999,fin);
	/* fprintf(stderr,"%s",s); */
	if(*s=='#'){
		/* comment - look for keyword */
		p=s;
		/* skip whitespace */
		while(1){
			p++;
			if(*p==' ')continue;
			if(*p=='\t')continue;
			break;
			}
		if(!strncasecmp(p,"GPS start:",10)){
			sscanf(p+10,"%Ld",gps);
			if(!check_intervals(segment_list, *gps)){
				fclose(fin);
				return -1;
				}
			if(check_intervals(veto_segment_list, *gps)>0){
				fclose(fin);
				return -1;
				}
			}
		}
	}
fgets(s,1999,fin);
if(strncmp(s,"EgiB",4)){
	fprintf(stderr,"File \"%s\" gps=%lld has wrong endianness.\n", filename,*gps);
	fclose(fin);
	return -1;
	}
fseek(fin, startbin*4,SEEK_CUR);
if(fread(data,4,count,fin)<count){
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	fclose(fin);
	return -1;
	}
fclose(fin);
return 0;
}

int get_geo_range(char *filename, long startbin, long count, float *data, INT64 *gps)
{
FILE *fin;
char s[PATH_MAX];
char *p;
REAL8 a, timebase;
INT4 b, bin_start, nbins;
REAL4 *tmp;
float factor;
long i;
errno=0;
fin=fopen(filename,"r");
if(errno==ENOENT)return -1; /* no such file */
if(fin==NULL){
	fprintf(stderr,"Error reading band %ld-%ld bins from file \"%s\":",startbin,startbin+count,filename);
	perror("");
	return -1;	
	}
/* read header */	
/* Key */	
fread(&a, sizeof(a), 1, fin);
if(a!=1.0){
	fprintf(stderr,"Cannot read file \"%s\": wrong endianness\n", filename);
	return -1;
	}
/* gps */
fread(&b, sizeof(b), 1, fin);
*gps=b;

if(!check_intervals(segment_list, *gps)){
	fclose(fin);
	return -1;
	}
if(check_intervals(veto_segment_list, *gps)>0){
	fclose(fin);
	return -1;
	}


/* skip nsec */
fread(&b, sizeof(b), 1, fin);
/* timebase */
fread(&timebase, sizeof(a), 1, fin);

fread(&bin_start, sizeof(bin_start), 1, fin);
fread(&nbins, sizeof(nbins), 1, fin);

tmp=do_alloc(count*2, sizeof(*tmp));

fseek(fin, (startbin-bin_start)*8,SEEK_CUR);
if(fread(tmp,4,count*2,fin)<count*2){
	fprintf(stderr,"Not enough data in file \"%s\" gps=%lld.\n",filename,*gps);
	free(tmp);
	fclose(fin);
	return -1;
	}
/* reverse normalization applied to geo format files */
if(timebase < 0) {
	factor=1.0; /* make_sft_op did not apply normalization .. */
	fprintf(stderr,"** Timebase is negative, assuming unnormalized data\n");
	fprintf(LOG,"** Timebase is negative, assuming unnormalized data\n");
	} else {
	factor=(0.5*1800.0*16384.0)/nbins; /* use fixed normalization for 1800 sec SFTs .. */
	}
factor*=factor; /* square it */
for(i=0;i<count;i++){
	data[i]=(tmp[2*i]*tmp[2*i]+tmp[2*i+1]*tmp[2*i+1])*factor;
	}
free(tmp);
fclose(fin);
return 0;
}


void read_directory(char *prefix, long first,long last, 
	long first_bin,long bin_count,
	long *nsegments, float **power, INT64 **gps)
{
char s[PATH_MAX];
long i;
int (*get_range)(char *filename, long startbin, long count, float *data, INT64 *gps);

fprintf(LOG,"power data size: %f MB\n",(last-first+1)*bin_count*sizeof(**power)/(1024.0*1024.0));
*power=do_alloc(last-first+1,bin_count*sizeof(**power));
*gps=do_alloc(last-first+1,sizeof(**gps));
*nsegments=0;

if(!strcasecmp("GEO", args_info.input_format_arg)){
	get_range=get_geo_range;
	} else
if(!strcasecmp("Power", args_info.input_format_arg)){
	get_range=get_power_range;
	} else {
	fprintf(stderr,"Unknown input format option: \"%s\"\n", args_info.input_format_arg);
	exit(-1);
	}

fprintf(stderr,"Reading files %s*:", prefix);
for(i=first;i<=last;i++){
	snprintf(s,PATH_MAX, args_info.input_munch_arg, prefix, i);
	fprintf(stderr," %ld",i);
	if(!get_range(s,first_bin,bin_count,*power+(*nsegments)*bin_count,(*gps)+(*nsegments))){
		fprintf(stderr,"(%Ld)",(*gps)[*nsegments]);
		(*nsegments)++;
		}
	}
fprintf(stderr,"\nRead %ld files\n", *nsegments);
/* free memory we did not use */
if(*nsegments<(last-first+1)){
	*power=realloc(*power,(*nsegments)*bin_count*sizeof(**power));
	*gps=realloc(*gps,(*nsegments)*sizeof(**gps));
	}
}


int clear_name_dat(char *name)
{
return !regexec(&write_dat, name, 0, NULL, 0);
}

int clear_name_png(char *name)
{
return !regexec(&write_png, name, 0, NULL, 0);
}

void dump_shorts(char *name, short *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "shorts: %s\n", name);
}

void dump_ints(char *name, int *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "ints: %s\n", name);
}

void dump_floats(char *name, float *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "floats: %s\n", name);
}

void dump_doubles(char *name, double *x, long count, long step)
{
FILE *fout;
long i;
char s[PATH_MAX];

if(!clear_name_dat(name))return;

snprintf(s,PATH_MAX,"%s%s", output_dir, name);
fout=fopen(s, "w");
if(fout==NULL){
	fprintf(FILE_LOG, "Could not open file \"%s\" for writing.\n",
		name);
	return;
	}
/* important special case */
if(step==1)fwrite(x, sizeof(*x), count,fout);
	else
	for(i=0;i<count;i++)fwrite(x+i, sizeof(*x), 1,fout);
fclose(fout);
fprintf(FILE_LOG, "doubles: %s\n", name);
}


void init_ephemeris(void)
{
LALStatus status={level:0, statusPtr:NULL};

memset(&ephemeris, 0, sizeof(ephemeris));

if(!strcasecmp("LHO", args_info.detector_arg)){
	detector=lalCachedDetectors[LALDetectorIndexLHODIFF];
	} else 
if(!strcasecmp("LLO", args_info.detector_arg)){
	detector=lalCachedDetectors[LALDetectorIndexLLODIFF];
	} else {
	fprintf(stderr,"Unrecognized detector site: \"%s\"\n", args_info.detector_arg);
	exit(-1);
	}
fprintf(LOG,"detector  : %s (%s)\n", args_info.detector_arg, detector.frDetector.name);
ephemeris.ephiles.earthEphemeris=earth_ephemeris;
ephemeris.ephiles.sunEphemeris=sun_ephemeris;
LALInitBarycenter(&status, &ephemeris);
TESTSTATUS(&status);

#if 0
detectorvel_inputs.detector=detector;
detectorvel_inputs.edat=&ephemeris;
#endif
fprintf(stderr,"Successfully initialized ephemeris data\n");
}

void get_AM_response(INT64 gps, float latitude, float longitude, float orientation,
	float *plus, float *cross)
{
LALStatus status={level:0, statusPtr:NULL};
INT4 tmp_leapsecs;
LALGPSandAcc gps_and_acc;
LALLeapSecFormatAndAcc leapsec_info={LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
LALSource source;
LALDetAndSource det_and_source={NULL, NULL};
LALDetAMResponse response;

memset(&gps_and_acc, 0, sizeof(gps_and_acc));
gps_and_acc.accuracy=leapsec_info.accuracy;
gps_and_acc.gps.gpsSeconds=gps; 
gps_and_acc.gps.gpsNanoSeconds=0;

memset(&source, 0, sizeof(source));
source.equatorialCoords.system=COORDINATESYSTEM_EQUATORIAL;
source.orientation=orientation;
source.equatorialCoords.longitude=longitude;
source.equatorialCoords.latitude=latitude;

det_and_source.pDetector=&detector;
det_and_source.pSource=&source;

/* TODO : find out why DetAMResponse does not care about leap seconds, but
   DetectorVel does */

LALComputeDetAMResponse(&status, &response, &det_and_source, &gps_and_acc);
TESTSTATUS(&status);

*cross=response.cross;
*plus=response.plus;
}

void get_detector_vel(INT64 gps, float *velocity)
{
LALStatus status={level:0, statusPtr:NULL};
REAL8 det_velocity[3];
INT4 tmp_leapsecs;
LALGPSandAcc gps_and_acc;
LALLeapSecFormatAndAcc  leapsec_info={LALLEAPSEC_GPSUTC, LALLEAPSEC_STRICT};
/* LIGOTimeGPS gps; */
int i;

memset(&gps_and_acc, 0, sizeof(gps_and_acc));
gps_and_acc.accuracy=leapsec_info.accuracy;
gps_and_acc.gps.gpsSeconds=gps; 
gps_and_acc.gps.gpsNanoSeconds=0;
 
LALLeapSecs(&status, &tmp_leapsecs, &(gps_and_acc.gps), &leapsec_info);
TESTSTATUS(&status);
ephemeris.leap = (INT2)tmp_leapsecs; 

LALDetectorVel(&status, det_velocity, &(gps_and_acc.gps), detector, &ephemeris);
TESTSTATUS(&status);

#if 0
fprintf(stderr,"powerflux: det_velocity=(%g,%g,%g)\n", 
	det_velocity[0],
	det_velocity[1],
	det_velocity[2]
	);
fprintf(stderr,"gps=%d (nano=%d)\n",gps_and_acc.gps.gpsSeconds, gps_and_acc.gps.gpsNanoSeconds);
fprintf(stderr,"detector=%s\n", detector.frDetector.name);
fprintf(stderr,"powerflux leap=%d nE=%d nS=%d dE=%g dS=%g\n", 
	ephemeris.leap,
	ephemeris.nentriesE,
	ephemeris.nentriesS,
	ephemeris.dtEtable,
	ephemeris.dtStable);
#endif

for(i=0;i<3;i++)velocity[i]=det_velocity[i];
}


/* there are count*GRID_FIT_COUNT coefficients */
void get_whole_sky_AM_response(INT64 *gps, long count, float orientation, float **coeffs_plus, float **coeffs_cross, long *size)
{
long i, j, k;
SKY_GRID *sample_grid=NULL;
float plus, cross;

gsl_multifit_linear_workspace *workspace=NULL;
gsl_matrix *X=NULL, *cov=NULL;
gsl_vector *y_plus=NULL, *y_cross=NULL, *c=NULL;
double chisq;

fprintf(stderr,"Computing whole sky AM response\n");
fprintf(LOG, "AM coeffs size: %f MB\n", 2*count*GRID_FIT_COUNT*sizeof(*coeffs_plus)/(1024.0*1024.0));
*size=count*GRID_FIT_COUNT;
*coeffs_plus=do_alloc(*size, sizeof(**coeffs_plus));
*coeffs_cross=do_alloc(*size, sizeof(**coeffs_cross));

sample_grid=make_rect_grid(8,5);

workspace=gsl_multifit_linear_alloc(sample_grid->npoints, GRID_FIT_COUNT);

y_plus=gsl_vector_alloc(sample_grid->npoints);
y_cross=gsl_vector_alloc(sample_grid->npoints);
c=gsl_vector_alloc(GRID_FIT_COUNT);

cov=gsl_matrix_alloc(GRID_FIT_COUNT, GRID_FIT_COUNT);
X=gsl_matrix_alloc(sample_grid->npoints, GRID_FIT_COUNT);

for(k=0;k<count;k++){
	for(i=0;i<sample_grid->npoints;i++){
		get_AM_response(gps[k]+900, 
			sample_grid->latitude[i], sample_grid->longitude[i], 
			orientation,
			&plus, &cross);
		gsl_vector_set(y_plus, i, plus);
		gsl_vector_set(y_cross, i, cross);
		
		for(j=0;j<GRID_FIT_COUNT;j++){
			gsl_matrix_set(X, i, j, sample_grid->e[j+GRID_FIT_START][i]);
			}
		}
	gsl_multifit_linear(X, y_plus, c, cov, &chisq, workspace);
	if(chisq>1e-12){
		fprintf(stderr,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, plus polarization - aborting !\n", 
			chisq, k, gps[k]);
		fprintf(LOG,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, plus polarization - aborting !\n", 
			chisq, k, gps[k]);
		exit(-1);
		}
	for(j=0;j<GRID_FIT_COUNT;j++){
		(*coeffs_plus)[k*GRID_FIT_COUNT+j]=gsl_vector_get(c, j);
		}

	gsl_multifit_linear(X, y_cross, c, cov, &chisq, workspace);
	if(chisq>1e-12){
		fprintf(stderr,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, cross polarization - aborting !\n", 
			chisq, k, gps[k]);
		fprintf(LOG,"** Sky grid approximation fault: non-zero chisq %g when computing gps[%d]=%lld, cross polarization - aborting !\n", 
			chisq, k, gps[k]);
		exit(-1);
		}
	for(j=0;j<GRID_FIT_COUNT;j++){
		(*coeffs_cross)[k*GRID_FIT_COUNT+j]=gsl_vector_get(c, j);
		}
	}

gsl_matrix_free(X);
gsl_matrix_free(cov);
gsl_vector_free(y_plus);
gsl_vector_free(y_cross);
gsl_vector_free(c);

gsl_multifit_linear_free(workspace);
free_grid(sample_grid);
}

/* there are count*GRID_FIT_COUNT coefficients */
void verify_whole_sky_AM_response(INT64 *gps, long count, float orientation,  SKY_GRID *sample_grid, float *coeffs_plus, char *name)
{
int i,j;
long offset;
float plus, cross;
float max_err,err;
gsl_rng *rng=NULL;
rng=gsl_rng_alloc(gsl_rng_default);

max_err=0;
for(i=0;i<count;i++){
	/* 20 points per segment ought to be enough */
	for(j=0;j<20;j++){
		/* test in random grid points */
		offset=floor(sample_grid->npoints*gsl_rng_uniform(rng));
		get_AM_response(gps[i]+900, 
			sample_grid->latitude[offset], sample_grid->longitude[offset], 
			orientation,
			&plus, &cross);
		err=fabs(plus*plus-AM_response(i, sample_grid, offset, coeffs_plus));
		if(err>max_err)max_err=err;
		}
	}
fprintf(stderr, "%s AM coeffs error: %g\n", name, max_err);
fprintf(LOG, "%s AM coeffs error: %g\n", name, max_err);
gsl_rng_free(rng);
}
