#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/SkyCoordinates.h>

REAL8TimeSeries * XLALSignalDetectorStrainREAL8TimeSeries( REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, SkyPosition *position, REAL8 psi, LALDetector *detector );

REAL8TimeSeries * XLALInjectionREAL8TimeSeries( REAL8TimeSeries *h, LIGOTimeGPS *start, REAL8 deltaT, UINT4 length, COMPLEX8FrequencySeries *response );
