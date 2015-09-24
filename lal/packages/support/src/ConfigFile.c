/*
 * Copyright (C) 2010 Larne Pekowsky
 * Copyright (C) 2004, 2005 Reinhard Prix
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/** \file
 * \ingroup ConfigFile
 * \author Reinhard Prix
 *
 * \brief General-purpose routines for config-file reading.
 */

/**
 * \addtogroup ConfigFile

   \par Description

This module provides routines for reading formatted
config-files containing definitions of the form <tt>variable = value</tt>.
The general syntax is somewhat similar to the one provided by the
perl-module <tt>ConfigParser</tt> (cf.
http://www.python.org/doc/current/lib/module-ConfigParser.html)

Comments are allowed using either '<tt>\#</tt>', '<tt>;</tt>' or <tt>%</tt>.
You can also use line-continuation  using a '<tt>\\</tt>' at the end of the line.
Also note that comment-signs '<tt>\#;%</tt>' within double-quotes '<tt>"..."</tt>'
are <em>not</em> treated as comment-characters.  The general syntax is best illustrated
using a simple example:
\code
# comment line
var1 = 1.0    ; you can also comment using semi-colons
somevar = some text.\
        You can also use\
        line-continuation
   var3 = 4      # whatever that means
note = "this is also possible, and # here does nothing"
a_switch = true	 #possible values: 0,1,true,false,yes,no, case insensitive
%% etc etc.
\endcode

Note that TABS generally get replaced by a single space, which can be
useful in the case of line-continuation (see example). All leading and
trailing spaces in are ignore (except within double-quotes).

The general approach of reading from such a config-file, is to first
call XLALParseDataFile() which loads and pre-parses the contents of the
config-file into the structure LALParsedDataFile. Then one can read in
config-variables either using one of the type-strict custom-wrappers
<tt>XLALReadConfig<TYPE>Variable()</tt> or the general-purpose reading function
XLALReadConfigVariable().


A boolean variable read by XLALReadConfigBOOLVariable() can have any of the values
<tt>{"1", "0", "yes", "no", "true", "false"}</tt>, where the comparison is done
<em>case-insensitively</em>, i.e. you can also use "True" or "FALSE"....


If one wishes a tight sytnax for the config-file, one can check
that there are no illegal entries in the config-file. This is done
by checking at the end that all config-file entries have been
successfully parsed, using:
XLALCheckConfigReadComplete(), where \a strictness is either
CONFIGFILE_WARN or CONFIGFILE_ERROR.
In the first case only a warning is issued, while in the second it is
treated as a LAL-error if some config-file entries have not been
read-in. (The use of this function is optional).


The configfile-data should be freed at the end using
XLALDestroyParsedDataFile().

\par Uses
\code
XLALCHARReadSequence()
XLALCreateTokenList()       XLALDestroyTokenList()
XLALCalloc()                XLALMalloc()             XLALFree()
XLALPrintError()            LALOpenDataFile()        fclose()
\endcode

\par Notes

XLALReadConfigSTRINGVariable() and XLALReadConfigSTRINGVariable() are not
the same as using <tt>"%s"</tt> as a format string, as they read the
<em>rest</em> of the logical line (excluding comments) as a string.


In the case of XLALReadConfigSTRINGVariable(), the required
memory is allocated and has to be freed by the caller, while for
XLALReadConfigSTRINGVariable() the caller has to provide a
CHARVector of length \f$N\f$, which defines the maximum length of
string to be read.


\note instead of using these functions directly, it might be
more convenient to use the UserInput.c infrastructure

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <errno.h>

#if HAVE_SYS_STAT_H
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

/* #include <ctype.h> */  /* don't use this, as it binds us to GLIBC_2.3 symbols!! */

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/StreamInput.h>
#include <lal/LogPrintf.h>

#include <lal/ConfigFile.h>

NRCSID( CONFIGFILEC, "$Id$");

extern INT4 lalDebugLevel;

#define FMT_STRING "string"    /* reading in quoted strings needs some special treatment */
#define WHITESPACE " \t"

#define TRUE   (1==1)
#define FALSE  (1==0)

/* FIXME: prototype for function in StreamSequenceInput.m4 */
int XLALCHARReadSequence( CHARSequence **sequence, FILE *stream );

/* local prototypes */
static void cleanConfig (CHARSequence *text);
CHAR my_tolower (CHAR in);
/* ctype replacements w/o locale */
static int TOLOWER(int c);
static const char upper_chars[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char lower_chars[] = "abcdefghijklmnopqrstuvwxyz";



/** Parse an ASCII data-file into a pre-cleaned array of lines.
 *
 * The cleaning gets rid of comments ('#', ';'), empty lines,
 * and performs line-continuation if '\' is found at EOL
 */
int
XLALParseDataFile (LALParsedDataFile **cfgdata, /**< [out] pre-parsed data-file lines */
		   const CHAR *fname)		/**< [in] name of config-file to be read */
{
  const char *fn = __func__;
  CHARSequence *rawdata = NULL;
  FILE *fp;
  int err = 0;  /* error code */
#if HAVE_STAT
  struct stat stat_out;
#endif

  if (*cfgdata != NULL) {
    XLALPrintError ("%s:" CONFIGFILEH_MSGENONULL, fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if (fname == NULL) {
    XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

#if HAVE_STAT
  if (  stat ( fname, &stat_out ) )
    {
      XLALPrintError ( "%s: Could not stat data-file: `%s` : \n\n", fn, fname, strerror(errno) );
      XLALPrintError( CONFIGFILEH_MSGEFILE);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if ( ! S_ISREG(stat_out.st_mode)  )
    {
      XLALPrintError ( "%s: '%s' does not seem to be a regular file!\n", fn, fname );
      XLALPrintError ( CONFIGFILEH_MSGEFILE );
      XLAL_ERROR ( fn, XLAL_EDOM );
    }
#endif

  if ( (fp = LALOpenDataFile (fname)) == NULL) {
    XLALPrintError ( "%s: Could not open data-file: `%s`\n\n", fn, fname);
    XLALPrintError (CONFIGFILEH_MSGEFILE);
    XLAL_ERROR ( fn, XLAL_ESYS );
  }

  err = XLALCHARReadSequence (&rawdata, fp);
  fclose (fp);
  if (err)
    return err;

  if (rawdata == NULL) {
    XLALPrintError( "%s:" CONFIGFILEH_MSGEFILE, fn);
    XLAL_ERROR ( fn, XLAL_ESYS );
  }

  /* get rid of comments and do line-continuation */
  cleanConfig (rawdata);

  if ( (*cfgdata = XLALCalloc (1, sizeof(LALParsedDataFile))) == NULL) {
    XLALPrintError ("%s:" CONFIGFILEH_MSGEMEM, fn );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

  /* parse this into individual lines */
  err = XLALCreateTokenList (&((*cfgdata)->lines), rawdata->data, "\n");
  XLALFree (rawdata->data);
  XLALFree (rawdata);

  if (err) {
    XLALPrintError ("%s: XLALCreateTokenList() failed.\n", fn );
    XLALFree (*cfgdata);
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  /* initialize the 'wasRead' flags for the lines */
  if ( (*cfgdata)->lines->nTokens )
    {
      if ( ((*cfgdata)->wasRead = XLALCalloc(1, (*cfgdata)->lines->nTokens * sizeof( (*cfgdata)->wasRead[0]))) == NULL )
	{
	  XLALFree ((*cfgdata)->lines);
	  XLALPrintError ( "%s:" CONFIGFILEH_MSGEMEM, fn );
	  XLAL_ERROR ( fn, XLAL_ENOMEM );
	}
    }
  else
    (*cfgdata)->wasRead = NULL;

  return XLAL_SUCCESS;

} /* XLALParseDataFile() */


/** Free memory associated with a LALParsedDataFile structure.
 */
int
XLALDestroyParsedDataFile (LALParsedDataFile **cfgdata)	/**< [in/out] config-file data */
{
  const char *fn = __func__;
  int err = 0;  /* error code */

  if ( ! cfgdata || ! *cfgdata || ! (*cfgdata)->lines ) {
    XLALPrintError( "%s:" CONFIGFILEH_MSGENULL, fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  err = XLALDestroyTokenList ( &((*cfgdata)->lines) );
  if (err) {
    XLALPrintError("%s: XLALDestroyTokenList() failed.\n", fn );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  if ( (*cfgdata)->wasRead )
    XLALFree ( (*cfgdata)->wasRead );

  XLALFree ( *cfgdata );

  *cfgdata = NULL;

  return XLAL_SUCCESS;

} /* XLALDestroyParsedDataFile() */


/** Function to determine whether a given section secName exists in the parsed
 * config-file contents cfgdata.
 *
 * \note: this function tolerates NULL input as secName, cfgdata, or cfgdata->lines,
 * in which case the answer is simply 'FALSE'.
 */
int
XLALConfigSectionExists ( const LALParsedDataFile *cfgdata,  	/**< [in] pre-parsed config-data */
                          const CHAR *secName)			/**< [in] section-name to read */
{
  const char *fn = __func__;
  UINT4 i;
  size_t sec_searchlen = 0;

  /* This traps coding errors in the calling routine. */
  /* If there's no config file, or no section, then   */
  /* the section isn;t in the config file, return 0   */
  if (secName == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn );
      return FALSE;
    }

  sec_searchlen = strlen(secName);

  if (cfgdata == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn );
      return FALSE;
    }

  if (cfgdata->lines == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn );
      return FALSE;
    }

  for (i = 0; i < cfgdata->lines->nTokens; i++)
    {
      /* Is this the start of a new section? */
      if (cfgdata->lines->tokens[i][0] == '[')
        {
          /* If we're looking for a particular section, is this it? */
          if (strncmp(cfgdata->lines->tokens[i] + 1, secName, sec_searchlen) == 0)
            {
              return TRUE;
            }
         }
    }

  return FALSE;

} /* XLALConfigSectionExists() */


/** Parser for config-file: can read config-variables of the form
 *	VARIABLE [=:] VALUE.
 * Input is a TokenList containing the 'logical' lines of the cleaned config-file
 *
 * - <tt>param->varName</tt> is the name of the config-variable to read
 * - <tt>param->fmt</tt>     is the format string to use for reading
 *
 * \note a special format-string is FMT_STRING, which means read the whole remaining line
 *   which is different from "%s"! (reads only one word)
 *   In this case, this also does the memory-allocation!
 *
 */
int
XLALReadConfigVariable ( void *varp,			  /**< [out] result gets written here! */
                         const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
                         const LALConfigVar *param,	  /**< [in]  var-name, fmt-string, strictness */
                         BOOLEAN *wasRead)		  /**< [out] did we succeed in reading? */
{
  const char *fn = __func__;
  CHAR *found    = NULL;
  INT2 ret = 0;

  UINT4 i;
  INT4 linefound     = -1;
  INT4 section_found = -1;

  size_t len;
  size_t searchlen = strlen (param->varName);
  size_t sec_searchlen = 0;

  if (param->secName == 0)
    {
      /* If we haven't been asked for a section then we want the
         "default" section, which starts at the top of the file */
      sec_searchlen = 0;
      section_found = 1;
    } else {
      sec_searchlen = strlen(param->secName);
    }

  /* This traps coding errors in the calling routine. */
  if (cfgdata == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (cfgdata->lines == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (cfgdata->wasRead == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (varp == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (param->varName == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (param->fmt == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  *wasRead = FALSE;

  /* let's look for the variable-name in the token-list (has to at beginning of line!) */
  for (i = 0; i < cfgdata->lines->nTokens; i++)
    {
      /* Is this the start of a new section? */
      if (cfgdata->lines->tokens[i][0] == '[')
        {
          /* If we're looking for a particular section, is this it? */
          if (sec_searchlen > 0 && strncmp(cfgdata->lines->tokens[i] + 1, param->secName, sec_searchlen) == 0)
            {
              section_found = i;
            }
          else
            section_found = -1;  /* We might have moved out of the right section */
         }
      else if (section_found > -1) /* Not section start, are we in the one we want? */
        {
          len = strcspn (cfgdata->lines->tokens[i], WHITESPACE "=:");	/* get length of variable-name */
          if (len == 0)
	        {			/* malformed token-list */
	          XLALPrintError ("%s:" CONFIGFILEH_MSGETOKENS, fn);
	          XLAL_ERROR ( fn, XLAL_EDOM );
	        }

          /* pre-select based on length of variable-name */
          if (len != searchlen)
            continue;

          /* same len, but are they identical ? */
          if (strncmp (param->varName, cfgdata->lines->tokens[i], len) == 0)
	        {
	          found = cfgdata->lines->tokens[i] + len;
	          found += strspn (found, WHITESPACE "=:");	/* skip all whitespace and define-chars */
	          linefound = i;
	          break;		/* ok, we've found it */
	        }
        }           /* if section found */
    }				/* for lines */

  if (!found)
    {
      switch (param->strictness)
	{
	case CONFIGFILE_IGNORE:
	  return 0;
	  break;
	case CONFIGFILE_WARN:
	  if (lalDebugLevel & LALWARNING)
            {
              if (sec_searchlen > 0)
                XLALPrintError ("%s: Warning: Config-file variable '%s' in section '%s' was not found!\n",
                                fn, param->varName, param->secName);
              else
                XLALPrintError ("%s: Warning: Config-file variable '%s' was not found!\n", fn, param->varName );
            }
	  return 0;
	  break;
	case CONFIGFILE_ERROR:
	default:
          if (sec_searchlen > 0)
            XLALPrintError ( "%s: Error: Config-file variable '%s' in section '%s' was not found!\n", fn, param->varName, param->secName);
          else
            XLALPrintError ( "%s: Error: Config-file variable '%s' was not found!\n", fn, param->varName);
	  XLALPrintError (CONFIGFILEH_MSGEVAR);
	  XLAL_ERROR ( fn, XLAL_EDOM );
	  break;
	} /* switch (strictness) */

    } /* if not found */

  /* now read the value into the variable */

  /* reading a quoted string needs some special treatment: */
  if (!strcmp (param->fmt, FMT_STRING))
    {
      /* NOTE: varp here is supposed to be a pointer to CHAR* !! */
      CHAR **cstr = (CHAR **) varp;

      if (*cstr != NULL)
        {
          XLALPrintError ("%s:" CONFIGFILEH_MSGENONULL, fn );
          XLAL_ERROR ( fn, XLAL_EINVAL );
        }

      (*cstr) = (CHAR *) XLALMalloc (strlen (found) + 1);
      strcpy ((*cstr), found);
      ret = 1;
    }
  else				/* but the default case is just sscanf... */
    ret = sscanf (found, param->fmt, varp);

  if ((ret == 0) || (ret == EOF))
    {
      XLALPrintError ("%s: ERROR: Config-file variable %s was not readable using the format %s\n\n", fn, param->varName, param->fmt);
      XLALPrintError(CONFIGFILEH_MSGEFMT);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  /* ok, we have successfully read in the config-variable: let's make a note of it */
  cfgdata->wasRead[linefound] = 1;

  *wasRead = TRUE;

  return XLAL_SUCCESS;

} /* XLALReadConfigVariable() */


/** Type-specialization of generic reading-function XLALReadConfigVariable() to BOOLEAN variables.
 */
int
XLALReadConfigBOOLVariable (BOOLEAN *varp,		   /**< [out] variable to store result */
			    const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
			    const CHAR *secName,		       /**< [in] section-name to read */
			    const CHAR *varName,		       /**< [in] variable-name to read */
			    BOOLEAN *wasRead)			       /**< [out] did we succeed in reading? */
{
  const char *fn = __func__;
  CHAR *tmp = NULL;
  INT2 ret = -1;		/* -1 means no legal value has been parsed */
  INT2 ret2;
  *wasRead = FALSE;

  /* first read the value as a string */
  if ( XLALReadConfigSTRINGVariable (&tmp, cfgdata, secName, varName, wasRead) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALReadConfigSTRINGVariable() failed.\n", fn);
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  if (*wasRead && tmp)		/* if we read anything at all... */
    {
      /* get rid of case ambiguities */
      ret2 = XLALLowerCaseString (tmp);

      if (ret2)
        return ret2;

      /* try to parse it as a bool */
      if (!strcmp (tmp, "yes") || !strcmp (tmp, "true") || !strcmp (tmp, "1"))
        ret = 1;
      else if (!strcmp (tmp, "no") || !strcmp (tmp, "false") || !strcmp (tmp, "0"))
        ret = 0;
      else
        {
          XLALPrintError ("%s: illegal bool-value `%s`\n", fn, tmp);
          XLALPrintError (CONFIGFILEH_MSGEBOOL);
          XLALFree (tmp);
          XLAL_ERROR ( fn, XLAL_EINVAL );
        }

      XLALFree (tmp);

      if (ret != -1)		/* only set value of something has been found */
        {
          *varp = (BOOLEAN) ret;
          *wasRead = TRUE;
        }
    } /* if wasRead && tmp */


  return XLAL_SUCCESS;

} /* XLALReadConfigBOOLVariable() */


/** Type-specialization of generic reading-function LALReadConfigVariable() to INT4 variables.
 */
int
XLALReadConfigINT4Variable (INT4 *varp,
			    const LALParsedDataFile *cfgdata,
                            const CHAR *secName,
			    const CHAR *varName,
                            BOOLEAN *wasRead)
{
  LALConfigVar param = { 0, 0, 0, 0 };

  param.varName    = varName;
  param.fmt        = "%" LAL_INT4_FORMAT;
  param.strictness = CONFIGFILE_IGNORE;
  param.secName    = secName;

  return XLALReadConfigVariable ((void *) varp, cfgdata, &param, wasRead);

} /* XLALReadConfigINT4Variable() */


/** Type-specialization of generic reading-function LALReadConfigVariable() to REAL8 variables.
 */
int
XLALReadConfigREAL8Variable (REAL8 *varp,
			    const LALParsedDataFile *cfgdata,
			    const CHAR *secName,
			    const CHAR *varName,
			    BOOLEAN *wasRead)
{
  LALConfigVar param = {0,0,0,0};

  param.secName = secName;
  param.varName = varName;
  param.fmt = "%" LAL_REAL8_FORMAT;
  param.strictness = CONFIGFILE_IGNORE;

  return XLALReadConfigVariable( (void*) varp, cfgdata, &param, wasRead);

} /* XLALReadConfigREAL8Variable() */


/** Type-specialization of generic reading-function XLALReadConfigVariable() to
 * STRING variables
 * \note this means the rest of the line, NOT "%s"! (but excluding comments of course),
 * \par Note2: if string is quoted by ", everything within quotes is read,
 *       and the quotes are removed here
 *
 */
int
XLALReadConfigSTRINGVariable (CHAR ** varp,		/**< [out] string, allocated here! */
			      const LALParsedDataFile * cfgdata, /**< [in] pre-parsed config-data */
			      const CHAR * secName,	/**< [in] section-name to be read */
			      const CHAR * varName,	/**< [in] variable-name to be read */
			      BOOLEAN * wasRead)	/**< [out] did we succeed in reading? */
{
  const char *fn = __func__;
  LALConfigVar param = { 0, 0, 0, 0 };
  CHAR *str = NULL;
  CHAR *ret = NULL;

  if (*varp != NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENONULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  param.secName = secName;
  param.varName = varName;
  param.fmt = FMT_STRING;
  param.strictness = CONFIGFILE_IGNORE;

  if ( XLALReadConfigVariable ((void *) &str, cfgdata, &param, wasRead) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALReadConfigVariable() failed\n", fn );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  if (*wasRead && (str != NULL))
    {
      INT2 numQuotes = 0;
      CHAR *ptr = str;
      /* count number of quotation marks */
      while ((ptr = strchr (ptr, '"')))
	{
	  numQuotes++;
	  ptr++;
	}			/* while quotes found */

      /* check balanced quotes (don't allow escaping for now) */
      if ((numQuotes != 0) && (numQuotes != 2))
	{
	  XLALPrintError ("%s:" CONFIGFILEH_MSGESTRING, fn );
	  XLAL_ERROR ( fn, XLAL_EDOM );
	}
      if (numQuotes == 2)
	{
	  /* allowed only at end and beginning */
	  if ((str[0] != '"') || (str[strlen (str) - 1] != '"'))
	    {
	      XLALPrintError ("%s:" CONFIGFILEH_MSGESTRING, fn);
	      XLAL_ERROR ( fn, XLAL_EDOM );
	    }
	  /* quotes ok, now remove them */
	  if ((ret = XLALMalloc (strlen (str) - 2 + 1)) == NULL)
	    {
	      XLALPrintError ("%s:" CONFIGFILEH_MSGEMEM, fn);
	      XLAL_ERROR ( fn, XLAL_ENOMEM );
	    }
	  str[strlen (str) - 1] = 0;
	  strcpy (ret, str + 1);
	  XLALFree (str);
	}			/* if 2 quotation marks */
      else
	ret = str;		/* no quotes, just return string */

      *varp = ret;

    }				/* if wasRead */
  else
    *varp = NULL;

  return XLAL_SUCCESS;

} /* XLALReadConfigSTRINGVariable() */


/** Type-specialization of generic reading-function XLALReadConfigVariable() to
 * reading of <em>fixed-length</em> strings.
 * Another variant of string-reading:similar to ReadConfigSTRINGVariable(), but
 * here a fixed-size CHAR-array is used as input, no memory is allocated by
 * the function.
 * \note you have to provide the length of your string-array as input in <tt>varp->length</tt>
 * (this is basically a wrapper for ReadConfigSTRINGVariable())
 *
 * \par Note 2: the behaviour is similar to strncpy, i.e. we silently clip the
 *       string to the right length, BUT we also 0-terminate it properly.
 *       No error or warning is generated when clipping occurs!
 *
 * \par Note 3: at return, the value <tt>varp->length</tt> is set to the length of the
 *        string copied
 *
 */
int
XLALReadConfigSTRINGNVariable (CHARVector *varp, 	/**< [out] must be allocated! */
			      const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
			      const CHAR *secName,	/**< [in] section-name */
			      const CHAR *varName,	/**< [in] variable-name */
			      BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  const char *fn = __func__;
  CHAR *tmp = NULL;

  /* This traps coding errors in the calling routine. */
  if (varp == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if (varp->data == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  if ( XLALReadConfigSTRINGVariable (&tmp, cfgdata, secName, varName, wasRead) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALReadConfigSTRINGVariable() failed.\n", fn );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  if (*wasRead && tmp)
    {
      strncpy (varp->data, tmp, varp->length - 1);
      varp->data[varp->length-1] = '\0';
      XLALFree (tmp);
      varp->length = strlen (varp->data);
      *wasRead = TRUE;
    }
  else
    *wasRead = FALSE;

  return XLAL_SUCCESS;

} /* XLALReadConfigSTRINGNVariable() */


/** Check if all lines of config-file have been successfully read in
 * and issue a warning or error (depending on strictness) if not.
 */
int
XLALCheckConfigReadComplete (const LALParsedDataFile *cfgdata, 	/**< [in] config-file data */
                             ConfigStrictness strict)  		/**< [in] what to do if unparsed lines */
{
  const char *fn = __func__;
  UINT4 i;

  if (cfgdata == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }
  if (cfgdata->lines == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn);
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }
  if (cfgdata->wasRead == NULL)
    {
      XLALPrintError ("%s:" CONFIGFILEH_MSGENULL, fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  for (i=0; i < cfgdata->lines->nTokens; i++)
    {
      /* Don't require section headers to be marked read */
      /* This has the effect of considering a file to have */
      /* been read completely if every value in every section */
      /* has been read. */
      if (cfgdata->lines->tokens[i][0] == '[')
        continue;

      if (cfgdata->wasRead[i] == 0)
	{
	  switch (strict)
	    {
	    case CONFIGFILE_IGNORE:
	      continue;
	    case CONFIGFILE_WARN:
	      XLALPrintError ( "%s: Warning: Ignoring unknown config-file entry '%s'.\n", fn, cfgdata->lines->tokens[i] );
	      continue;
	    case CONFIGFILE_ERROR:
	    default:
	      XLALPrintError ( "%s: ERROR: config-file entry #%d has not been read!\n", fn, i);
	      XLALPrintError ( "Line was: '%s'\n", cfgdata->lines->tokens[i]);
              XLALPrintError (CONFIGFILEH_MSGEUNKNOWN);
              XLAL_ERROR ( fn, XLAL_EDOM );
	    } /* switch strict */
	} /* if some line not read */

    } /* for i < lines */

  return XLAL_SUCCESS;

} /* XLALCheckConfigReadComplete() */


/** Helper function:  turn a string into lowercase without using locale-functions.
 */
int
XLALLowerCaseString (CHAR *string)	/**< [in/out] string to convert */
{
  const char *fn = __func__;
  UINT4 i;

  if (string == NULL)
    {
      XLALPrintError ( "%s:" CONFIGFILEH_MSGENULL, fn );
      XLAL_ERROR ( fn, XLAL_EINVAL );
    }

  for (i=0; i < strlen (string); i++)
    string[i] = TOLOWER( string[i] );

  return XLAL_SUCCESS;

} /* XLALLowerCaseString() */

/*----------------------------------------------------------------------
 * tolower() replacement w/o locale
 *----------------------------------------------------------------------*/
static int
TOLOWER(int c)
{
  if (c) {
    char *p = strchr(upper_chars, c);

    if (p) {
      c = lower_chars[p - upper_chars];
    }
  }
  return c;
} /* TOLOWER() */

/*----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------
 * cleanConfig(): do some preprocessing on the config-file, namely 'erase'
 * all comments by '\n', and glue '\'-continued lines
 *----------------------------------------------------------------------*/
void
cleanConfig (CHARSequence *text)
{
  size_t len;
  CHAR *ptr, *ptr2, *eol;
  BOOLEAN inQuotes = 0;

  if ( !text )
    return;

  /*----------------------------------------------------------------------
   * RUN 1: clean out comments, by replacing them by '\n'
   */
  ptr = text->data;
  while ( *ptr )
    {
      if ( (*ptr) == '\"' )
	inQuotes = !inQuotes;	/* flip state */

      if ( ((*ptr) == '#') || ( (*ptr) == ';') || ( (*ptr) == '%') )
	if ( !inQuotes )	/* only consider as comments if not quoted */
	  {
	    len = strcspn (ptr, "\n");
	    memset ( (void*)ptr, '\n', len);
	  }

      ptr ++;

    } /* while *ptr */

  /*----------------------------------------------------------------------
   * RUN 2: do line-gluing when '\' is found at end-of-line
   */
  ptr = text->data;
  while ( (ptr = strchr(ptr, '\\')) != NULL )
    {
      if ( ptr[1] == '\n' )
	{
	  /* ok, now it gets a bit tricky: to avoid getting spurious spaces from
	   * the line-continuation, we shift the rest of the file forward by 2 positions
	   * to nicely fit to the previous line...
	   */
	  len = strlen (ptr+2);
	  memmove(ptr, ptr+2, len+1);	/* move the whole rest (add +1 for '\0') */
	}
    } /* while '\' found in text */

  /*----------------------------------------------------------------------
   * RUN 3: turn all tabs into single spaces..
   */
  ptr = text->data;
  while ( (ptr = strchr(ptr, '\t')) != NULL )
    *ptr = ' ';

  /*----------------------------------------------------------------------
   * RUN 4: get rid of initial and trailing whitespace (replace it by '\n')
   */
  ptr = text->data;
  while (ptr < (text->data + text->length -1) )
    {
      eol = strchr (ptr, '\n'); /* point to end-of-line */

      len = strspn (ptr, WHITESPACE);
      if (len) memset ( (void*)ptr, '\n', len);

      if (eol != NULL)
	ptr = eol;
      else
	ptr = strchr (ptr, '\0'); /* or end of file */

      /* clean away all trailing whitespace of last line*/
      ptr2 = ptr - 1;
      while ( *ptr2 == ' ' )
	*ptr2-- = '\n';

      /* step to next line */
      ptr += 1;
    }

  return;

} /* cleanConfig() */


/* ========== DEPRECATED LAL INTERFACE FUNCTIONS, which have been replaced by XLAL functions,
 * These functions are just wrappers around the XLAL functions
 */

/** Deprecated LAL interface wrapper to XLALParseDataFile()
 *
 */
void
LALParseDataFile (LALStatus *status,		/**< pointer to LALStatus structure */
		  LALParsedDataFile **cfgdata, 	/**< [out] pre-parsed data-file lines */
		  const CHAR *fname)		/**< [in] name of config-file to be read */
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALParseDataFile (cfgdata, fname) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALParseDataFile() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALLoadConfigFile() */


/** Deprecated LAL interface wrapper to XLALDestroyParsedDataFile()
 */
void
LALDestroyParsedDataFile (LALStatus *status,		/**< pointer to LALStatus structure */
			  LALParsedDataFile **cfgdata)	/**< [in/out] config-file data */
{
  const char *fn = __func__;
  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALDestroyParsedDataFile (cfgdata) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALDestroyParsedDataFile() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALDestroyConfigData() */




/** Deprecated LAL interface wrapper to XLALReadConfigVariable()
 */
void
LALReadConfigVariable (LALStatus *status,		/**< pointer to LALStatus structure */
		       void *varp, 			/**< [out] result gets written here! */
		       const LALParsedDataFile *cfgdata,/**< [in] pre-parsed config-data */
		       const LALConfigVar *param,	/**< [in]  var-name, fmt-string, strictness */
		       BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALReadConfigVariable ( varp,	cfgdata, param, wasRead ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALReadConfigVariable() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALReadConfigVariable() */


/** Deprecated LAL interface wrapper to XLALReadConfigBOOLVariable()
 */
void
LALReadConfigBOOLVariable (LALStatus *status,		/**< pointer to LALStatus structure */
			   BOOLEAN *varp, 		 /**< [out] variable to store result */
			   const LALParsedDataFile *cfgdata,/**< [in] pre-parsed config-data */
			   const CHAR *varName,		 /**< [in] variable-name to read */
			   BOOLEAN *wasRead)		 /**< [out] did we succeed in reading? */
{
  const char *fn = __func__;
  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALReadConfigBOOLVariable(varp, cfgdata, NULL, varName, wasRead ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALReadConfigBOOLVariable() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALReadConfigBOOLVariable() */


/** Deprecated LAL interface wrapper to XLALReadConfigINT4Variable()
 */
void
LALReadConfigINT4Variable (LALStatus *status,
			   INT4 *varp,
			   const LALParsedDataFile *cfgdata,
			   const CHAR *varName,
			   BOOLEAN *wasRead)
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALReadConfigINT4Variable ( varp, cfgdata, NULL, varName, wasRead ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALReadConfigINT4Variable() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALReadConfigINT4Variable() */

/** Deprecated LAL interface wrapper to XLALReadConfigREAL8Variable()
 */
void
LALReadConfigREAL8Variable (LALStatus *status,
			    REAL8 *varp,
			    const LALParsedDataFile *cfgdata,
			    const CHAR *varName,
			    BOOLEAN *wasRead)
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALReadConfigREAL8Variable ( varp, cfgdata, NULL, varName, wasRead ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALReadConfigREAL8Variable() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALReadConfigREAL8Variable() */

/** Deprecated LAL interface wrapper to XLALReadConfigSTRINGVariable()
 */
void
LALReadConfigSTRINGVariable (LALStatus *status,		/**< pointer to LALStatus structure */
			     CHAR **varp, 		/**< [out] string, allocated here! */
			     const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
			     const CHAR *varName,	/**< [in] variable-name to be read */
			     BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALReadConfigSTRINGVariable ( varp, cfgdata, NULL, varName, wasRead ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALReadConfigSTRINGVariable() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALReadConfigSTRINGVariable() */


/** Deprecated LAL interface wrapper to XLALReadConfigSTRINGNVariable()
 */
void
LALReadConfigSTRINGNVariable (LALStatus *status,	/**< pointer to LALStatus structure */
			      CHARVector *varp, 	/**< [out] must be allocated! */
			      const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
			      const CHAR *varName,	/**< [in] variable-name */
			      BOOLEAN *wasRead)		/**< [out] did we succeed in reading? */
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALReadConfigSTRINGNVariable(varp, cfgdata, NULL, varName, wasRead ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALReadConfigSTRINGNVariable() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALReadConfigSTRINGNVariable() */

/** Deprecated LAL interface wrapper to XLALCheckConfigReadComplete()
 */
void
LALCheckConfigReadComplete (LALStatus *status,			/**< pointer to LALStatus structure */
			    const LALParsedDataFile *cfgdata, 	/**< [in] config-file data */
			    ConfigStrictness strict)  		/**< [in] what to do if unparsed lines */
{
  const char *fn = __func__;

  INITSTATUS( status, fn, CONFIGFILEC );

  if ( XLALCheckConfigReadComplete ( cfgdata, strict ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: call to XLALCheckConfigReadComplete() failed with code %d\n", fn, xlalErrno );
    ABORT ( status, CONFIGFILEH_EXLAL, CONFIGFILEH_MSGEXLAL );
  }

  RETURN (status);

} /* LALCheckConfigReadComplete() */