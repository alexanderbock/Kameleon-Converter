/******************************************************************************
* Copyright 1996-2013 United States Government as represented by the
* Administrator of the National Aeronautics and Space Administration.
* All Rights Reserved.
******************************************************************************/
/****************************************************************************** 
*
*  SPDF/CDF                       CDF_TIME_TT2000 utility routines for C.
*
*  Version 1.0, 15-Mar-11, ADNET systems
*
*  Modification history:
*
*   V1.0  15-Mar-11, M Liu      Initial version (for CDF V3.3.2)
*   V2.0  15-Feb-12, M Liu      Initial version (for CDF V3.4.0)
******************************************************************************/

#include "cdflib.h"
#include "cdflib64.h"

#define YearWithin(a)		((a >= 1708) && (a <= 2291))
#define JulianDateJ2000_12h     2451545
#define J2000Since0AD12h        730485 /* 2000-1-1 since 0-1-1 */
#define J2000Since0AD12hSec     63113904000.0
#define J2000Since0AD12hMilsec  63113904000000.0
                                /* 86400000.0 * 730485 */
#define J2000LeapSeconds        32.0
#define dT                      32.184
#define dTinNanoSecs            32184000000LL
#define MJDbase                 2400000.5
#define SECinNanoSecs           1000000000LL
#define SECinNanoSecsD          1000000000.0
#define DAYinNanoSecs           86400000000000LL
#define HOURinNanoSecs          3600000000000LL
#define MINUTEinNanoSecs        60000000000LL
#define T12hinNanoSecs          43200000000000LL
/* Julian days for 1707-09-22 and 2292-04-11, the valid TT2000 range. */
#define JDY17070922             2344793
#define JDY22920411             2558297

/* Number of Delta(dAT) expressions before leap seconds were introduced */
#define NERA1 14

static int tableChecked = 0;
static int openCDF64s = 0;
static int toPlus = 0;
static long currentDay = -1;
static double currentLeapSeconds;
static double currentJDay;
static int fromFile = 0;

/* Dates, Delta(AT)s and drift rates */
static double  LTS [][6]= {
/*      year month day    delta       drift      drift     */
      { 1960,  1,   1,  1.4178180,   37300.0, 0.0012960 },
      { 1961,  1,   1,  1.4228180,   37300.0, 0.0012960 },
      { 1961,  8,   1,  1.3728180,   37300.0, 0.0012960 },
      { 1962,  1,   1,  1.8458580,   37665.0, 0.0011232 },
      { 1963, 11,   1,  1.9458580,   37665.0, 0.0011232 },
      { 1964,  1,   1,  3.2401300,   38761.0, 0.0012960 },
      { 1964,  4,   1,  3.3401300,   38761.0, 0.0012960 },
      { 1964,  9,   1,  3.4401300,   38761.0, 0.0012960 },
      { 1965,  1,   1,  3.5401300,   38761.0, 0.0012960 },
      { 1965,  3,   1,  3.6401300,   38761.0, 0.0012960 },
      { 1965,  7,   1,  3.7401300,   38761.0, 0.0012960 },
      { 1965,  9,   1,  3.8401300,   38761.0, 0.0012960 },
      { 1966,  1,   1,  4.3131700,   39126.0, 0.0025920 },
      { 1968,  2,   1,  4.2131700,   39126.0, 0.0025920 },
      { 1972,  1,   1, 10.0,             0.0, 0.0       },
      { 1972,  7,   1, 11.0,             0.0, 0.0       },
      { 1973,  1,   1, 12.0,             0.0, 0.0       },
      { 1974,  1,   1, 13.0,             0.0, 0.0       },
      { 1975,  1,   1, 14.0,             0.0, 0.0       },
      { 1976,  1,   1, 15.0,             0.0, 0.0       },
      { 1977,  1,   1, 16.0,             0.0, 0.0       },
      { 1978,  1,   1, 17.0,             0.0, 0.0       },
      { 1979,  1,   1, 18.0,             0.0, 0.0       },
      { 1980,  1,   1, 19.0,             0.0, 0.0       },
      { 1981,  7,   1, 20.0,             0.0, 0.0       },
      { 1982,  7,   1, 21.0,             0.0, 0.0       },
      { 1983,  7,   1, 22.0,             0.0, 0.0       },
      { 1985,  7,   1, 23.0,             0.0, 0.0       },
      { 1988,  1,   1, 24.0,             0.0, 0.0       },
      { 1990,  1,   1, 25.0,             0.0, 0.0       },
      { 1991,  1,   1, 26.0,             0.0, 0.0       },
      { 1992,  7,   1, 27.0,             0.0, 0.0       },
      { 1993,  7,   1, 28.0,             0.0, 0.0       },
      { 1994,  7,   1, 29.0,             0.0, 0.0       },
      { 1996,  1,   1, 30.0,             0.0, 0.0       },
      { 1997,  7,   1, 31.0,             0.0, 0.0       },
      { 1999,  1,   1, 32.0,             0.0, 0.0       },
      { 2006,  1,   1, 33.0,             0.0, 0.0       },
      { 2009,  1,   1, 34.0,             0.0, 0.0       },
      { 2012,  7,   1, 35.0,             0.0, 0.0       }
     };

/* Number of Delta(AT) in static table */
static const int NDAT = sizeof (LTS) / sizeof (LTS[0]);
static double **LTD = NULL;
static long long *NST = NULL;
static int entryCnt;

static double JulianDay12h (long, long, long);
static void DatefromJulianDay (double, long *, long *, long *);
static int ValidateYMD (long, long , long);
static char *MonthToken (long);
static int MonthNumber (char *);
static int ScanUTCstring (char *);
static double LeapSecondsfromYMD (long, long, long);
static double LeapSecondsfromJ2000 (long long, int *);
static void LoadLeapSecondsTable ();
static void LoadLeapNanoSecondsTable ();
static void EPOCHbreakdownTT2000 (double, long *, long *, long *, long *,
                                  long *, long *);

static char *MonthToken (long month)
{
  switch (month) {
    case 1: return "Jan";
    case 2: return "Feb";
    case 3: return "Mar";
    case 4: return "Apr";
    case 5: return "May";
    case 6: return "Jun";
    case 7: return "Jul";
    case 8: return "Aug";
    case 9: return "Sep";
    case 10: return "Oct";
    case 11: return "Nov";
    case 12: return "Dec";
  }
  return "???";
}

static int MonthNumber (char *month)
{
  if (StrStrIgCaseX (month, "jan")) return 1;
  else if (StrStrIgCaseX (month, "feb")) return 2;
  else if (StrStrIgCaseX (month, "mar")) return 3;
  else if (StrStrIgCaseX (month, "apr")) return 4;
  else if (StrStrIgCaseX (month, "may")) return 5;
  else if (StrStrIgCaseX (month, "jun")) return 6;
  else if (StrStrIgCaseX (month, "jul")) return 7;
  else if (StrStrIgCaseX (month, "aug")) return 8;
  else if (StrStrIgCaseX (month, "sep")) return 9;
  else if (StrStrIgCaseX (month, "oct")) return 10;
  else if (StrStrIgCaseX (month, "nov")) return 11;
  else if (StrStrIgCaseX (month, "dec")) return 12;
  else return 0;
}


static double JulianDay12h (long y, long m, long d)
{
  if (m == 0) m = 1;
  return (double) (367*y-7*(y+(m+9)/12)/4-3*((y+(m-9)/7)/100+1)/4+275*m/9+d+1721029);
}

static void DatefromJulianDay (double julday, long *y, long *m, long *d)
{
  long i, j, k, l, n;
  l=(long) (julday+68569);
  n=4*l/146097;
  l=l-(146097*n+3)/4;
  i=4000*(l+1)/1461001;
  l=l-1461*i/4+31;
  j=80*l/2447;
  k=l-2447*j/80;
  l=j/11;
  j=j+2-12*l;
  i=100*(n-49)+i+l;
  *y = i;
  *m = j;
  *d = k;
  return; 
}

static int ScanUTCstring (char *string) 
{
  int len;
  len = (int) strlen (string);
  if ((len == TT2000_3_STRING_LEN) ||
      (len > 20 && string[10] == 'T')) return 3;
  else if ((len == TT2000_0_STRING_LEN) || 
           (len > 20 && string[11] == ' ')) return 0; 
  else if ((len == TT2000_1_STRING_LEN) ||
           (len > 9 && string[8] == '.')) return 1;
  else if ((len == TT2000_2_STRING_LEN) || (len > 7))  return 2;
  return -1;
}

static int ValidateYMD (long yy, long mm, long dd)
{
  double jday;
  if (yy <= 0 || mm < 0 || dd < 0) return 0;
  /* Y-M-D should be in the 1707-09-22 and 2292-04-11 range. */
  jday = JulianDay12h(yy, mm, dd);
  if (jday < JDY17070922 || jday > JDY22920411) return 0;
  return 1; 
}

void EPOCHbreakdownTT2000 (double epoch, long *year, long *month, long *day,
                           long *hour, long *minute, long *second)
{
  long jd,i,j,k,l,n;
  double second_AD, minute_AD, hour_AD, day_AD;
  second_AD = epoch;
  minute_AD = second_AD / 60.0;
  hour_AD = minute_AD / 60.0;
  day_AD = hour_AD / 24.0;

  jd = (long) (1721060 + day_AD);
  l=jd+68569;
  n=4*l/146097;
  l=l-(146097*n+3)/4;
  i=4000*(l+1)/1461001;
  l=l-1461*i/4+31;
  j=80*l/2447;
  k=l-2447*j/80;
  l=j/11;
  j=j+2-12*l;
  i=100*(n-49)+i+l;

  *year = i;
  *month = j;
  *day = k;

  *hour   = (long) fmod (hour_AD, (double) 24.0);
  *minute = (long) fmod (minute_AD, (double) 60.0);
  *second = (long) fmod (second_AD, (double) 60.0);

  return;
}

static void LoadLeapSecondsTable ()
{
  if (!tableChecked) {
    char *table;
    FILE *leaptable = NULL;
    int im, ix;
    tableChecked = 1;
    table = CDFgetLeapSecondsTableEnvVar();
    if (table != NULL) {
      leaptable = fopen (table, "r");
      if (leaptable != NULL) {
        int togo = 1;
        char line[81];
        int count = 0;
        long yy, mm, dd;
        while (fgets(line, 80, leaptable)) {
          if (line[0] == ';') continue;
          ++count;
        }
        rewind(leaptable);
        LTD = (double **) cdf_AllocateMemory (count * sizeof (double *), NULL);
        ix = 0;
        while (fgets(line, 80, leaptable)) {
          if (line[0] == ';') continue;
          LTD[ix] = (double *) cdf_AllocateMemory(6 * sizeof(double), NULL);
          im = sscanf (line, "%ld %ld %ld %lf %lf %lf", &yy, &mm, &dd,
                       &(LTD[ix][3]), &(LTD[ix][4]), &(LTD[ix][5]));
          if (im != 6) {
            int iz;
            for (iz = 0; iz < ix; ++iz)
              cdf_FreeMemory (LTD[iz], NULL);
            cdf_FreeMemory (LTD, NULL);
            togo = 0;
            break;
          }
          LTD[ix][0] = (double) yy;
          LTD[ix][1] = (double) mm;
          LTD[ix][2] = (double) dd;
          ++ix;
        }
        fclose (leaptable);
        if (togo == 1) {
          entryCnt = count;
          fromFile = 1;
        }
      }
    }
    if (fromFile == 0) {
      LTD = (double **) cdf_AllocateMemory (NDAT * sizeof (double *), NULL);
      for (ix = 0; ix < NDAT; ++ix) {
        LTD[ix] = (double *) cdf_AllocateMemory(6 * sizeof(double), NULL);
        LTD[ix][0] = LTS[ix][0];
        LTD[ix][1] = LTS[ix][1];
        LTD[ix][2] = LTS[ix][2];
        LTD[ix][3] = LTS[ix][3];
        LTD[ix][4] = LTS[ix][4];
        LTD[ix][5] = LTS[ix][5];
      }
      entryCnt = NDAT;
    }
  }
}

static void LoadLeapNanoSecondsTable () 
{
  int ix;
  if (LTD == NULL) LoadLeapSecondsTable ();
  NST = (long long *) cdf_AllocateMemory (entryCnt * sizeof (long long), NULL);
  for (ix = NERA1; ix < entryCnt; ++ix)
     NST[ix] = CDF_TT2000_from_UTC_parts(LTD[ix][0], LTD[ix][1],
                                         LTD[ix][2], 0.0, 0.0, 0.0,
                                         0.0, 0.0, 0.0);
}

static double LeapSecondsfromYMD (long iy, long im, long id)
{
  int i, j;
  long m, n;
  double da;
  if (LTD == NULL) LoadLeapSecondsTable();
  j = -1;
  m = 12 * iy + im;
  for (i = entryCnt-1; i >=0; --i) {
    n = (long) (12 * LTD[i][0] + LTD[i][1]);
    if (m >= n) {
      j = i;
      break;
    }
  }
  if (j == -1) return 0.0;
  da = LTD[j][3];
  /* If pre-1972, adjust for drift. */
  if (j < NERA1) {
    double jda;
    jda = JulianDay12h (iy, im, id);
    da += ((jda - MJDbase) - LTD[j][4]) * LTD[j][5];
  }
  return da;
}

static double LeapSecondsfromJ2000 (long long nanosecs, int *leapSecond)
{
  int i, j;
  double da;
  *leapSecond = 0;
  if (NST == NULL) LoadLeapNanoSecondsTable();
  j = -1;
  for (i = entryCnt-1; i >=NERA1; --i) {
    if (nanosecs >= NST[i]) {
      j = i;
      if (i < (entryCnt - 1)) {
        /* Check for time following on leap second (second = 60). */
        if ((nanosecs + 1000000000L) >= NST[i+1]) {
          *leapSecond = 1;
        }
      }
      break;
    }
  }
  if (j == -1) return 0.0; /* Pre 1972 .... do it later.... */
  da = LTD[j][3];
  return da;
}

/******************************************************************************
* CDFgetLeapSecondsTableEnvVar.
******************************************************************************/
    
VISIBLE_PREFIX char *CDFgetLeapSecondsTableEnvVar ()
{
  char *table = NULL;
#if defined(vms)
    table = getenv("CDF$LEAPSECONDSTABLE");
#else
    table = getenv("CDF_LEAPSECONDSTABLE");
#endif
    return table;
}

/******************************************************************************
* CDFClearLeapSecondsTable.
******************************************************************************/
    
VISIBLE_PREFIX void CDFClearLeapSecondsTable ()
{
  int ix;
  if (LTD != NULL) {
    for (ix = 0; ix < entryCnt; ++ix)
       cdf_FreeMemory (LTD[ix], NULL);
    cdf_FreeMemory (LTD, NULL);
    if (NST != NULL) cdf_FreeMemory (NST, NULL);
    LTD = NULL;
    NST = NULL;
  }
  tableChecked  = 0;
}

/******************************************************************************
* CDFfillLeapSecondsTable.
******************************************************************************/

VISIBLE_PREFIX void CDFfillLeapSecondsTable (double **table)
{
  int ix;
  if (LTD == NULL) LoadLeapSecondsTable();
  for (ix = 0; ix < NDAT; ++ix) {
    table[ix][0] = LTD[ix][0];
    table[ix][1] = LTD[ix][1];
    table[ix][2] = LTD[ix][2];
    table[ix][3] = LTD[ix][3];
    table[ix][4] = LTD[ix][4];
    table[ix][5] = LTD[ix][5];
  }
}

/******************************************************************************
* CDFgetRowsinLeapSecondsTable.
******************************************************************************/

VISIBLE_PREFIX int CDFgetRowsinLeapSecondsTable ()
{
  if (LTD == NULL) LoadLeapSecondsTable();
  return entryCnt;
}

/******************************************************************************
* CDFgetLastDateinLeapSecondsTable.
******************************************************************************/
#if defined(vms)      
VISIBLE_PREFIX void CDFgetLastDateinLeapSecondsTBL (long *year, long *month, 
                                                    long *day)
#else
VISIBLE_PREFIX void CDFgetLastDateinLeapSecondsTable (long *year, long *month,
                                                      long *day)
#endif
{   
  if (LTD == NULL) LoadLeapSecondsTable();
  *year = (long) LTD[entryCnt-1][0];
  *month = (long) LTD[entryCnt-1][1];
  *day = (long) LTD[entryCnt-1][2];
}     
      
/******************************************************************************
* CDFgetLeapSecondsTableStatus.
******************************************************************************/
      
VISIBLE_PREFIX int CDFgetLeapSecondsTableStatus ()
{   
  if (LTD == NULL) LoadLeapSecondsTable();
  return fromFile;
}
      
/******************************************************************************
* CDF_TT2000_to_UTC_parts.
******************************************************************************/

#if defined(STDARG)
VISIBLE_PREFIX void CDF_TT2000_to_UTC_parts (long long nanoSecSinceJ2000,
                                             double *ly, double *lm,
                                             double *ld, ...)
#else
VISIBLE_PREFIX void CDF_TT2000_to_UTC_parts (long long nanoSecSinceJ2000,
                                             double *ly, double *lm,
                                             double *ld, va_alist)
va_dcl
#endif
{
  double epoch, *tmp, tmp1, dat0;
  long long t2, t3, secSinceJ2000;
  long nansec, milsec;
  double *ye, *ym, *da, *opt[6];
  long ye1, mo1, da1, ho1, mi1, se1, ml1, ma1, na1;
  int  ix, leapSec;
  va_list Va;
  ye = ly;
  ym = lm;
  da = ld;
#if defined(STDARG)
  va_start (Va, ld);
#else
  VA_START (Va);
#endif
  ix = 0;
  tmp = va_arg (Va, double *);
  while (tmp != TT2000NULL) {
    opt[ix] = tmp;
    ++ix;
    if (ix == 6) break;
    tmp = va_arg (Va, double *);
  }
  va_end (Va);
  if (nanoSecSinceJ2000 == FILLED_TT2000_VALUE) {
    *ly = 9999.0;
    *lm = 12.0;
    *ld = 31.0;
    if (ix > 0) *opt[0] = 23.0;
    if (ix > 1) *opt[1] = 59.0;
    if (ix > 2) *opt[2] = 59.0;
    if (ix > 3) *opt[3] = 999.0;
    if (ix > 4) *opt[4] = 999.0;
    if (ix > 5) *opt[5] = 999.0;
    return;
/*
  } else if (nanoSecSinceJ2000 == ILLEGAL_TT2000_VALUE) {
    *ly = 1710.0;
    *lm = 1.0;
    *ld = 1.0;
    if (ix > 0) *opt[0] = 0.0;
    if (ix > 1) *opt[1] = 0.0;
    if (ix > 2) *opt[2] = 0.0;
    if (ix > 3) *opt[3] = 0.0;
    if (ix > 4) *opt[4] = 0.0;
    if (ix > 5) *opt[5] = 0.0;
    return;
*/
  } else if (nanoSecSinceJ2000 == DEFAULT_TT2000_PADVALUE) {
    *ly = 0.0;
    *lm = 1.0;
    *ld = 1.0;
    if (ix > 0) *opt[0] = 0.0;
    if (ix > 1) *opt[1] = 0.0;
    if (ix > 2) *opt[2] = 0.0;
    if (ix > 3) *opt[3] = 0.0;
    if (ix > 4) *opt[4] = 0.0;
    if (ix > 5) *opt[5] = 0.0;
    return;
  }
  toPlus = 0;
  t3 = nanoSecSinceJ2000;
  dat0 = LeapSecondsfromJ2000 (nanoSecSinceJ2000, &leapSec);
/*
  nanoSecSinceJ2000 -= dTinNanoSecs;
  nanoSecSinceJ2000 += T12hinNanoSecs;
  t2 = nanoSecSinceJ2000;
  secSinceJ2000 = (long long) ((double)nanoSecSinceJ2000/SECinNanoSecsD);
  nansec = (long) (nanoSecSinceJ2000 - secSinceJ2000 * SECinNanoSecs);
  if (nansec < 0) {
    nansec = SECinNanoSecs + nansec;
    --secSinceJ2000;
  }
*/
  if (nanoSecSinceJ2000 > 0) { /* try to avoid overflow (substraction first) */
    secSinceJ2000 = (long) ((double)nanoSecSinceJ2000/SECinNanoSecsD);
    nansec = (long) (nanoSecSinceJ2000 - secSinceJ2000 * SECinNanoSecs);
    secSinceJ2000 -= 32; /* secs portion in dT */
    secSinceJ2000 += 43200; /* secs in 12h */
    nansec -= 184000000L; /* nanosecs portion in dT */
  } else { /* try to avoid underflow (addition first) */
    nanoSecSinceJ2000 += T12hinNanoSecs; /* 12h in nanosecs */
    nanoSecSinceJ2000 -= dTinNanoSecs /* dT in nanosecs */;
    secSinceJ2000 = (long) ((double)nanoSecSinceJ2000/SECinNanoSecsD);
    nansec = (long) (nanoSecSinceJ2000 - secSinceJ2000 * SECinNanoSecs);
  }
  if (nansec < 0) {
    nansec = SECinNanoSecs + nansec;
    --secSinceJ2000;
  }
  t2 = secSinceJ2000 * SECinNanoSecs + nansec;
  if (dat0 > 0.0) { /* Post-1972.... */
    secSinceJ2000 -= (long long) dat0;
    epoch = (double) J2000Since0AD12hSec + secSinceJ2000;
    if (leapSec == 0) 
      EPOCHbreakdownTT2000 (epoch, &ye1, &mo1, &da1, &ho1, &mi1, &se1);
    else {
      /* second is at 60.... bearkdown function can't handle 60 so make it
         59 first and then add 1 second back. */
      epoch -= 1.0;
      EPOCHbreakdownTT2000 (epoch, &ye1, &mo1, &da1, &ho1, &mi1, &se1);
      se1 += 1;
    }
  } else { /* Pre-1972.... */
    long long tmpNanosecs;
    epoch = (double) secSinceJ2000 + J2000Since0AD12hSec;
    /* First guess */
    EPOCHbreakdownTT2000 (epoch, &ye1, &mo1, &da1, &ho1, &mi1, &se1);
    tmpNanosecs = CDF_TT2000_from_UTC_parts ((double) ye1, (double) mo1,
                                             (double) da1, (double) ho1,
                                             (double) mi1, (double) se1,
                                             0.0, 0.0, (double) nansec);
    if (tmpNanosecs != t3) {
      long long tmpx, tmpy;
      dat0 = LeapSecondsfromYMD (ye1, mo1, da1);
      tmpx = t2 - (long long) (dat0 * SECinNanoSecs);
      tmpy = (long long) ((double)tmpx/SECinNanoSecsD);
      nansec = (long) (tmpx - tmpy * SECinNanoSecs);
      if (nansec < 0) {
        nansec = SECinNanoSecs + nansec;
        --tmpy;
      }
      epoch = (double) tmpy + J2000Since0AD12hSec;
      /* Second guess */
      EPOCHbreakdownTT2000 (epoch, &ye1, &mo1, &da1, &ho1, &mi1, &se1);
      tmpNanosecs = CDF_TT2000_from_UTC_parts ((double) ye1, (double) mo1,
                                               (double) da1, (double) ho1,
                                               (double) mi1, (double) se1,
                                               0.0, 0.0, (double) nansec);
      if (tmpNanosecs != t3) {
        dat0 = LeapSecondsfromYMD (ye1, mo1, da1);
        tmpx = t2 - (long long) (dat0 * SECinNanoSecs);
        tmpy = (long long) ((double)tmpx/SECinNanoSecsD);
        nansec = (long) (tmpx - tmpy * SECinNanoSecs);
        if (nansec < 0) {
          nansec = SECinNanoSecs + nansec;
          --tmpy;
        }
        epoch = (double) tmpy + J2000Since0AD12hSec;
        /* One more determination */
        EPOCHbreakdownTT2000 (epoch, &ye1, &mo1, &da1, &ho1, &mi1, &se1);
      }
    }
  }
  if (se1 == 60) toPlus = 1; 
  ml1 = (long) (nansec / 1000000);
  tmp1 = nansec - 1000000 * ml1;
  if (ml1 > 1000) {
    ml1 -= 1000;
    se1 += 1;
  }
  ma1 = (long) (tmp1 / 1000);
  na1 = (long)  (tmp1 - 1000 * ma1);
  *ly = (double) ye1;
  *lm = (double) mo1;
  if (ix == 6) {
    *ld = (double) da1;
    *(opt[0]) = (double) ho1;
    *(opt[1]) = (double) mi1;
    *(opt[2]) = (double) se1;
    *(opt[3]) = (double) ml1;
    *(opt[4]) = (double) ma1;
    *(opt[5]) = (double) na1;
  } else if (ix == 5) {
    *ld = (double) da1;
    *(opt[0]) = (double) ho1;
    *(opt[1]) = (double) mi1;
    *(opt[2]) = (double) se1;
    *(opt[3]) = (double) ml1;
    *(opt[4]) = (double) ma1 + (na1 / 1000.0);
  } else if (ix == 4) {
    *ld = (double) da1;
    *(opt[0]) = (double) ho1;
    *(opt[1]) = (double) mi1;   
    *(opt[2]) = (double) se1;   
    *(opt[3]) = (double) ml1 + (ma1*1000.0 + na1) / 1000000.0;
  } else if (ix == 3) {
    *ld = (double) da1;
    *(opt[0]) = (double) ho1;
    *(opt[1]) = (double) mi1;
    tmp1 = ml1*1000000.0 + ma1*1000.0 + na1;
    *(opt[2]) = (double) se1 + tmp1 / 1000000000.0;
  } else if (ix == 2) {
    *ld = (double) da1;
    *(opt[0]) = (double) ho1;
    tmp1 = se1*1000000000.0 + ml1*1000000.0 + ma1*1000.0 + na1;
    *(opt[1]) = (double) mi1 + tmp1 / (60000000000.0+1000000000*toPlus);
  } else if (ix == 1) {
    *ld = (double) da1;
    tmp1 = mi1*60000000000.0 + se1*1000000000.0 + ml1*1000000.0 +
           ma1*1000.0 + na1;
    *(opt[0]) = (double) ho1 + tmp1 / (3600000000000.0+1000000000*toPlus);
  } else if (ix == 0) {
    tmp1 = ho1*3600000000000.0 + mi1*60000000000.0 + se1*1000000000.0 +
           ml1*1000000.0 + ma1*1000.0 + na1;
    *ld = (double) da1 + tmp1 / (86400000000000.0 + 1000000000*toPlus);
  }
}

/******************************************************************************
* CDF_TT2000_from_UTC_parts.
******************************************************************************/

#if defined(STDARG)
VISIBLE_PREFIX long long CDF_TT2000_from_UTC_parts (double yy, double mm, 
                                                    double dd, ...)
#else
VISIBLE_PREFIX long long CDF_TT2000_from_UTC_parts (double yy, double mm,
                                                    double dd, va_alist)
va_dcl
#endif
{
  double jd;
  long long subDayinNanoSecs, nanoSecSinceJ2000;
  long long t2, iy;
  double tmp, opt[6];
  va_list Va;
  int ix;
  double ly, lm, ld, lh, ln, ls, ll, lu, la;
  long lyl, lml, ldl, lhl, lnl, lsl, lll, lul, lal;
  double tmpy;
  long xx, xy, xm, xd;
#if defined(STDARG)
  va_start (Va, dd);
#else
  VA_START (Va);
#endif
  ix = 0;
  tmp = va_arg (Va, double);
  while (tmp != TT2000END) {
    opt[ix] = tmp;
    ++ix;
    if (ix == 6) break;
    tmp = va_arg (Va, double);
  }
  va_end (Va);
  ly = floor(yy);
  lm = floor(mm);
  ld = floor(dd);
  if (ix == 6) {
    if (opt[0] < 0.0 || opt[1] < 0.0 || opt[2] < 0.0 || opt[3] < 0.0 ||
        opt[4] < 0.0 || opt[5] < 0.0) return ILLEGAL_TT2000_VALUE;
    lh = opt[0];
    ln = opt[1];
    ls = opt[2];
    ll = opt[3];
    lu = opt[4];
    if ((dd-ld) != 0.0 || (lh-floor(lh)) != 0.0 || (ln-floor(ln)) != 0.0 ||
        (ls-floor(ls)) != 0.0 || (ll-floor(ll)) != 0.0 ||
        (lu-floor(lu)) != 0.0) return ILLEGAL_TT2000_VALUE;
    la = opt[5];
  } else if (ix == 5) {
    if (opt[0] < 0.0 || opt[1] < 0.0 || opt[2] < 0.0 || opt[3] < 0.0 ||
        opt[4] < 0.0) return ILLEGAL_TT2000_VALUE;
    lh = opt[0];
    ln = opt[1];
    ls = opt[2];
    ll = opt[3];
    if ((dd-ld) != 0.0 || (lh-floor(lh)) != 0.0 || (ln-floor(ln)) != 0.0 ||
        (ls-floor(ls)) != 0.0 || (ll-floor(ll)) != 0.0)
      return ILLEGAL_TT2000_VALUE;
    lu = floor(opt[4]);
    la = (opt[4] - lu) * 1000.0;
  } else if (ix == 4) {
    if (opt[0] < 0.0 || opt[1] < 0.0 || opt[2] < 0.0 || opt[3] < 0.0)
      return ILLEGAL_TT2000_VALUE;
    lh = opt[0];
    ln = opt[1];
    ls = opt[2];
    if ((dd-ld) != 0.0 || (lh-floor(lh)) != 0.0 || (ln-floor(ln)) != 0.0 ||
        (ls-floor(ls)) != 0.0) return ILLEGAL_TT2000_VALUE;
    ll = floor(opt[3]); 
    tmp = (opt[3] - ll) * 1000.0;
    lu = floor(tmp);
    la = (tmp - lu) * 1000.0;
  } else if (ix == 3) {
    if (opt[0] < 0.0 || opt[1] < 0.0 || opt[2] < 0.0)
      return ILLEGAL_TT2000_VALUE;
    lh = opt[0];
    ln = opt[1];
    if ((dd-ld) != 0.0 || (lh-floor(lh)) != 0.0 || (ln-floor(ln)) != 0.0)
      return ILLEGAL_TT2000_VALUE;
    ls = floor(opt[2]);
    tmp = (opt[2] - ls) * 1000.0;
    ll = floor(tmp); 
    tmp = (tmp - ll) * 1000.0;
    lu = floor(tmp);
    la = (tmp - lu) * 1000.0;
  } else if (ix == 2) {
    if (opt[0] < 0.0 || opt[1] < 0.0) return ILLEGAL_TT2000_VALUE;
    lh = opt[0];
    if ((dd-ld) != 0.0 || (lh-floor(lh)) != 0.0) return ILLEGAL_TT2000_VALUE;
    ln = floor(opt[1]);
    tmp = opt[1] - ln;
    if (tmp > 0.0) {
      tmp *= 60.0;
      ls = floor(tmp);
      tmp = (tmp - ls) * 1000.0;
      ll = floor(tmp);
      tmp = (tmp - ll) * 1000.0;
      lu = floor(tmp);
      la = (tmp - lu) * 1000.0;
    } else {
      ls = ll = lu = la = 0.0;
    }
  } else if (ix == 1) {
    if ((dd-ld) != 0.0) return ILLEGAL_TT2000_VALUE;
    if (opt[0] < 0.0) return ILLEGAL_TT2000_VALUE;
    tmp = opt[0];
    if (tmp > 0.0) {
      tmp *= 24.0;
      lh = floor(tmp);
      tmp = (tmp - lh) * 60.0;
      ln = floor(tmp);
      tmp = (tmp - ln) * 60.0;
      ls = floor(tmp);
      tmp = (tmp - ls) * 1000.0;
      ll = floor(tmp);
      tmp = (tmp - ll) * 1000.0;
      lu = floor(tmp);
      la = (tmp - lu) * 1000.0;
    } else {
      lh = ln = ls = ll = lu = la = 0.0;
    }
  } else { /* ix == 0 */
    tmp = dd - ld;
    if (tmp > 0.0) {
      tmp *= 24.0;
      lh = floor(tmp);
      tmp = (tmp - lh) * 60.0;
      ln = floor(tmp);
      tmp = (tmp - ln) * 60.0;
      ls = floor(tmp);
      tmp = (tmp - ls) * 1000.0;
      ll = floor(tmp);
      tmp = (tmp - ll) * 1000.0;
      lu = floor(tmp);
      la = (tmp - lu) * 1000.0;
    } else {
      lh = ln = ls = ll = lu = la = 0.0;
    }
  }
  lyl = lml = -999;
  if (la >= 1000.0) {
    double ad, ah, am, as, al, au;
    ad = floor(la / 86400000000000.0);
    la = la - ad * 86400000000000.0;
    ah = floor(la /3600000000000.0);
    la = la - ah * 3600000000000.0;
    am = floor(la / 60000000000.0);
    la = la - am * 60000000000.0;
    as = floor(la / 1000000000.0);
    la = la - as * 1000000000.0;
    al = floor(la /1000000.0); 
    la = la - al * 1000000.0;
    au = floor(la /1000.0);
    la = la - au * 1000.0;
    ld += ad;
    lh += ah;
    ln += am; 
    ls += as;
    ll += al; 
    lu += au;
    tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
    DatefromJulianDay (tmp, &lyl, &lml, &ldl);
  } 
  if (lu >= 1000.0) {
    double ad, ah, am, as, al;
    ad = floor(lu / 86400000000.0);
    lu = lu - ad * 86400000000.0;
    ah = floor(lu /3600000000.0);
    lu = lu - ah * 3600000000.0;
    am = floor(lu / 60000000);
    lu = lu - am * 60000000;
    as = floor(lu / 1000000);
    lu = lu - as * 1000000;
    al = floor(lu /1000);
    lu = lu - al * 1000;
    ld += ad;
    lh += ah;
    ln += am;
    ls += as;
    ll += al;
    tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
    DatefromJulianDay (tmp, &lyl, &lml, &ldl);
  }
  if (ll >= 1000.0) {
    double ad, ah, am, as;
    ad = floor(ll / 86400000);
    ll = ll - ad * 86400000;
    ah = floor(ll /3600000);
    ll = ll - ah * 3600000;
    am = floor(ll / 60000);
    ll = ll - am * 60000;
    as = floor(ll / 1000);
    ll = ll - as * 1000;
    ld += ad;
    lh += ah;
    ln += am;
    ls += as;
    tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
    DatefromJulianDay (tmp, &lyl, &lml, &ldl);
  }
  if (ls >= 60.0) {
    tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
    DatefromJulianDay (tmp+1, &xy, &xm, &xd);
    toPlus = LeapSecondsfromYMD(xy,xm,xd) - LeapSecondsfromYMD(ly,lm,ld);
    toPlus = floor(toPlus);
    if (ls >= (60.0+toPlus)) {
      double ad, ah, am;
      ad = floor(ls / (86400+toPlus));
      ls = ls - ad * (86400+toPlus);
      ah = floor(ls /(3600+toPlus));
      ls = ls - ah * (3600+toPlus);
      am = floor(ls / (60+toPlus));
      ls = ls - am * (60+toPlus);
      ld += ad;
      lh += ah;
      ln += am;
      tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
      DatefromJulianDay (tmp, &lyl, &lml, &ldl);
    }
  }
  if (ln >= 60.0) {
    double ad, ah;
    ad = floor(ln / 1440);
    ln = ln - (ad * 1440);
    ah = floor(ln / 60);
    ln = ln - ah * 60;
    ld += ad;
    lh += ah;
    tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
    DatefromJulianDay (tmp, &lyl, &lml, &ldl);
  }
  if (lh >= 24.0) {
    double ad;
    ad = floor(lh / 24.0); 
    lh = lh - ad * 24.0;
    ld += ad;
    tmp = JulianDay12h ((long) ly, (long) lm, (long) ld);
    DatefromJulianDay (tmp, &lyl, &lml, &ldl);
  }
  if (lyl == -999 && lml == -999) {
    lyl = (long) ly;
    lml = (long) lm;
    ldl = (long) ld;
  }
  lhl = (long) lh;
  lnl = (long) ln;
  lsl = (long) ls;
  lll = (long) ll;
  lul = (long) lu;
  lal = (long) la;
  if (lyl == 9999 && lml == 12 && ldl == 31 && lhl == 23 && lnl == 59 &&
      lsl == 59 && lll == 999 && lul == 999 && lal == 999)
    return FILLED_TT2000_VALUE;
  else if (lyl == 0 && lml == 1 && ldl == 1 && lhl == 0 && lnl == 0 &&
           lsl == 0 && lll == 0 && lul == 0 && lal == 0)
    return DEFAULT_TT2000_PADVALUE;
  if (!YearWithin(lyl) && !ValidateYMD(lyl,lml,ldl))
    return ILLEGAL_TT2000_VALUE;
  iy = 10000000 * lml + 10000 * ldl + lyl;
  if (iy != currentDay) {
    currentDay = iy;
    currentLeapSeconds = LeapSecondsfromYMD(lyl, lml, ldl);
    currentJDay = JulianDay12h(lyl,lml,ldl);
  }
  jd = currentJDay;
  jd = jd - JulianDateJ2000_12h;
  subDayinNanoSecs = lhl * HOURinNanoSecs + lnl * MINUTEinNanoSecs +
                     lsl * SECinNanoSecs + lll * 1000000 + lul * 1000 + lal;
  nanoSecSinceJ2000 = (long long) jd * DAYinNanoSecs + subDayinNanoSecs;
  t2 = (long long) (currentLeapSeconds * SECinNanoSecs);
  if (nanoSecSinceJ2000 < 0) {
    nanoSecSinceJ2000 += t2;
    nanoSecSinceJ2000 += dTinNanoSecs;
    nanoSecSinceJ2000 -= T12hinNanoSecs;
  } else {
    nanoSecSinceJ2000 -= T12hinNanoSecs;
    nanoSecSinceJ2000 += t2;
    nanoSecSinceJ2000 += dTinNanoSecs;
  }
  return nanoSecSinceJ2000;
}

/******************************************************************************
* CDF_TT2000_to_UTC_EPOCH.
******************************************************************************/

VISIBLE_PREFIX double CDF_TT2000_to_UTC_EPOCH (long long nanoSecSinceJ2000)
{
  double yy, mm, dd, hh, nn, ss, ll, uu, aa;
  double epoch;
  if (nanoSecSinceJ2000 == FILLED_TT2000_VALUE) return -1.0E31;
  else if (nanoSecSinceJ2000 == DEFAULT_TT2000_PADVALUE ||
           nanoSecSinceJ2000 == ILLEGAL_TT2000_VALUE) return 0.0;
  CDF_TT2000_to_UTC_parts (nanoSecSinceJ2000, &yy, &mm, &dd, &hh, &nn, &ss,
                           &ll, &uu, &aa);
  epoch = computeEPOCH ((long)yy, (long)mm, (long)dd, (long)hh, (long)nn,
                        (long)ss, (long)ll);
  return epoch;
}

/******************************************************************************
* CDF_TT2000_from_UTC_EPOCH.
******************************************************************************/

VISIBLE_PREFIX long long CDF_TT2000_from_UTC_EPOCH (double epoch)
{
  long long nanosecSinceJ2000;
  double millisec;
  long yy, mm, dd, hh, nn, ss, ll;
  if (epoch == -1.0E31 || epoch == -1.0E-31)
    return FILLED_TT2000_VALUE; 
  if (epoch == 0.0 || NegativeZeroReal8 (&epoch))
    return DEFAULT_TT2000_PADVALUE;
  EPOCHbreakdown (epoch, &yy, &mm, &dd, &hh, &nn, &ss, &ll);
  if (!YearWithin(yy) && !ValidateYMD(yy,mm,dd)) return ILLEGAL_TT2000_VALUE;
  return CDF_TT2000_from_UTC_parts ((double)yy, (double)mm, (double)dd,
                                    (double)hh, (double)nn, (double)ss,
                                    (double)ll, 0.0, 0.0);
}

/******************************************************************************
* CDF_TT2000_to_UTC_EPOCH16.
******************************************************************************/

VISIBLE_PREFIX double CDF_TT2000_to_UTC_EPOCH16 (long long nanoSecSinceJ2000,
                                                 double *epoch)
{
  double yy, mm, dd, hh, nn, ss, ll, uu, aa;
  double tmp;
  if (nanoSecSinceJ2000 == FILLED_TT2000_VALUE) {
    *epoch = -1.0E31;
    *(epoch+1) = -1.0E31;
    return 0.0;
  }
  else if (nanoSecSinceJ2000 == DEFAULT_TT2000_PADVALUE  ||
           nanoSecSinceJ2000 == ILLEGAL_TT2000_VALUE) {
    *epoch = 0.0;
    *(epoch+1) = 0.0;
    return 0.0;
  }
  CDF_TT2000_to_UTC_parts (nanoSecSinceJ2000, &yy, &mm, &dd, &hh, &nn, &ss,
                           &ll, &uu, &aa);
  tmp = computeEPOCH16 ((long)yy, (long)mm, (long)dd, (long)hh, (long)nn,
                        (long)ss, (long)ll, (long)uu, (long)aa, 0, epoch);
  return tmp;
}

/******************************************************************************
* CDF_TT2000_from_UTC_EPOCH16.
******************************************************************************/

VISIBLE_PREFIX long long CDF_TT2000_from_UTC_EPOCH16 (double *epoch)
{
  long nasec;
  long yy, mm, dd, hh, nn, ss, ll, uu, aa, pp;
  if (epoch[0] == 0.0 && epoch[1] == 0.0) return DEFAULT_TT2000_PADVALUE;
  if (NegativeZeroReal8 (&(epoch[0])) && NegativeZeroReal8 (&(epoch[1])))
    return DEFAULT_TT2000_PADVALUE;
  if (epoch[0] == -1.0E31 && epoch[1] == -1.0E31) return FILLED_TT2000_VALUE;
  if (epoch[0] == -1.0E-31 && epoch[1] == -1.0E-31) return FILLED_TT2000_VALUE;
  EPOCH16breakdown (epoch, &yy, &mm, &dd, &hh, &nn, &ss, &ll, &uu, &aa, &pp);
  if (!YearWithin(yy) && !ValidateYMD(yy,mm,dd)) return ILLEGAL_TT2000_VALUE;
  return CDF_TT2000_from_UTC_parts ((double)yy, (double)mm, (double)dd,
                                    (double)hh, (double)nn, (double)ss,
                                    (double)ll, (double)uu, (double)aa);
}

/******************************************************************************
* CDF_TT2000_to_UTC_string.
******************************************************************************/

#if defined(STDARG)
VISIBLE_PREFIX void CDF_TT2000_to_UTC_string (long long nanoSecSinceJ2000,
                                              char *string, ...)
#else
VISIBLE_PREFIX void CDF_TT2000_to_UTC_string (long long nanoSecSinceJ2000,
                                              char *string, va_alist)
va_dcl
#endif
{
  long ly, lm, ld, lh, ln, ls, ll, lu, la;
  double lyd, lmd, ldd, lhd, lnd, lsd, lld, lud, lad;
  int format;
  va_list Va;
#if defined(STDARG)
  va_start (Va, string);
#else
  VA_START (Va);
#endif
  format = va_arg (Va, int);
  if (format < 0 || format > 3) format = 3;
  va_end (Va);
  if (nanoSecSinceJ2000 == FILLED_TT2000_VALUE) {
    if (format == 0) {
      strncpy(string, "31-DEC-9999 23:59:59.999999999", 30);
      string[30] = '\0';
    } else if (format == 1) {
      strncpy(string, "99991231.9999999999", 19);
      string[19] = '\0';
    } else if (format == 2) {
      strncpy(string, "99991231235959", 14);
      string[14] = '\0';
    } else {
      strncpy(string, "9999-12-31T23:59:59.999999999", 29);
      string[29] = '\0';
    }
    return;
  } else if (nanoSecSinceJ2000 == DEFAULT_TT2000_PADVALUE) {
    if (format == 0) {
      strncpy(string, "01-JAN-0000 00:00:00.000000000", 30);
      string[30] = '\0';
    } else if (format == 1) {
      strncpy(string, "00000101.0000000000", 19);
      string[19] = '\0';
    } else if (format == 2) {
      strncpy(string, "00000101000000", 14);
      string[14] = '\0';
    } else {
      strncpy(string, "0000-01-01T00:00:00.000000000", 29);
      string[29] = '\0';
    }
    return;
  }
  CDF_TT2000_to_UTC_parts(nanoSecSinceJ2000, &lyd, &lmd, &ldd, &lhd, &lnd,
                          &lsd, &lld, &lud, &lad);
  ly = (long) lyd;
  lm = (long) lmd;
  ld = (long) ldd;
  lh = (long) lhd;
  ln = (long) lnd;
  ls = (long) lsd;
  ll = (long) lld;
  lu = (long) lud;
  la = (long) lad;
  if (format == 0) {
    snprintf(string, (size_t) TT2000_0_STRING_LEN+1, 
             "%2.2ld-%s-%4.4ld %2.2ld:%2.2ld:%2.2ld.%3.3ld%3.3ld%3.3ld",
             ld,MonthToken(lm),ly,lh,ln,ls,ll,lu,la);
  } else if (format == 1) {
    long milsecs, nansecs;
    double subday;
    long long subdayll;
    milsecs = 3600000 * lh + 60000 * ln + 1000 * ls + ll;
    nansecs = 1000 * lu + la;
    subday = (1000000.0 * (double) milsecs + nansecs)
             / (86400.0 * SECinNanoSecsD);
    subdayll = (long long) (subday * pow (10.0, 10));
    snprintf(string, (size_t) TT2000_1_STRING_LEN+1,
             "%4.4ld%2.2ld%2.2ld.%10.10lld", ly,lm,ld,subdayll);
  } else if (format == 2) {
    snprintf(string, (size_t) TT2000_2_STRING_LEN+1,
             "%4.4ld%2.2ld%2.2ld%2.2ld%2.2ld%2.2ld",ly,lm,ld,lh,ln,ls);
  } else if (format == 3) {
    snprintf(string, (size_t) TT2000_3_STRING_LEN+1,
             "%4.4ld-%2.2ld-%2.2ldT%2.2ld:%2.2ld:%2.2ld.%3.3ld%3.3ld%3.3ld",
             ly,lm,ld,lh,ln,ls,ll,lu,la);
  }
}

/******************************************************************************
* CDF_TT2000_from_UTC_string.
******************************************************************************/

VISIBLE_PREFIX long long CDF_TT2000_from_UTC_string (char *string)
{
  double fraction;
  long ly, lm, ld, lh, ln, ls, ll, lu, la;
  int len, format;

  while (string[0] == ' ') memmove (string, &string[1], strlen(string));
  format = ScanUTCstring(string);
  if (format == 0) {
    char moString[4];
    long tmp;
    if (sscanf(string,"%2ld-%c%c%c-%4ld %2ld:%2ld:%2ld.%9ld",
        &ld, &(moString[0]), &(moString[1]), &(moString[2]), &ly,
        &lh, &ln, &ls, &tmp) != 9) return ILLEGAL_TT2000_VALUE;
    moString[3] = '\0';
    lm = (long) MonthNumber(moString);
    if (ly == 9999 && lm == 12 && ld == 31 && lh == 23 && ln == 59 &&
        ls == 59 && tmp == 999999999)
      return FILLED_TT2000_VALUE;
    else if (ly == 0 && lm == 1 && ld == 1 && lh == 0 && ln == 0 &&
             ls == 0 && tmp == 0)
      return DEFAULT_TT2000_PADVALUE;
    if (tmp == 0) {
      ll = lu = la = 0;
    } else {
      char *dot = strrchr (string, '.');
      if (dot == NULL) {
        ll = lu = la = 0;
      } else {
        len = (int) strlen (dot + 1);
        if (len < 9) tmp = tmp * pow (10.0, 9-len);
        ll = (long) (tmp / 1000000);
        tmp = tmp - ll * 1000000;
        lu = (long) (tmp / 1000);
        la = (long) (tmp - lu * 1000);
      }
    }
    if (!YearWithin(ly) && !ValidateYMD(ly,lm,ld)) return ILLEGAL_TT2000_VALUE;
    return CDF_TT2000_from_UTC_parts ((double) ly, (double) lm, (double) ld,
                                      (double) lh, (double) ln, (double) ls,
                                      (double) ll, (double) lu, (double) la);
  } else if (format == 1) {
    long long tmp;
    if (sscanf(string,"%4ld%2ld%2ld.%lld",
        &ly, &lm, &ld, &tmp) != 4) return ILLEGAL_TT2000_VALUE;
    if (tmp == 0) 
      fraction = 0.0;
    else {
      char *dot = strrchr (string, '.');
      len = (int) strlen (dot + 1);
      fraction = tmp / pow(10.0, len);
    }
    if (ly == 9999 && lm == 12 && ld == 31 && tmp == 9999999999LL)
      return FILLED_TT2000_VALUE;
    else if (ly == 0 && lm == 1 && ld == 1 && tmp == 0)
      return DEFAULT_TT2000_PADVALUE;
    if (!YearWithin(ly) && !ValidateYMD(ly,lm,ld)) return ILLEGAL_TT2000_VALUE;
    return CDF_TT2000_from_UTC_parts ((double) ly, (double) lm, (double) ld,
                                      fraction, TT2000END);
  } else if (format == 2) {
    if (sscanf(string,"%4ld%2ld%2ld%2ld%2ld%2ld",
        &ly, &lm, &ld, &lh, &ln, &ls) != 6) return ILLEGAL_TT2000_VALUE;
    if (ly == 9999 && lm == 12 && ld == 31 && lh == 23 && ln == 59 &&
        ls == 59) return FILLED_TT2000_VALUE;
    else if (ly == 0 && lm == 1 && ld == 1 && lh == 0 && ln == 0 &&
             ls == 0) return DEFAULT_TT2000_PADVALUE;
    if (!YearWithin(ly) && !ValidateYMD(ly,lm,ld)) return ILLEGAL_TT2000_VALUE;
    return CDF_TT2000_from_UTC_parts ((double) ly, (double) lm, (double) ld,
                                      (double) lh, (double) ln, (double) ls,
                                      0.0, 0.0, 0.0);
  } else if (format == 3) {
    long tmp;
    if (sscanf(string,"%4ld-%2ld-%2ldT%2ld:%2ld:%2ld.%9ld",
        &ly, &lm, &ld, &lh, &ln, &ls, &tmp) != 7) return ILLEGAL_TT2000_VALUE;
    if (ly == 9999 && lm == 12 && ld == 31 && lh == 23 && ln == 59 &&
        ls == 59 && tmp == 999999999)
      return FILLED_TT2000_VALUE;
    else if (ly == 0 && lm == 1 && ld == 1 && lh == 0 && ln == 0 &&
             ls == 0 && tmp == 0)
      return DEFAULT_TT2000_PADVALUE;
    if (tmp == 0) {
      ll = lu = la = 0;
    } else {
      char *dot = strrchr (string, '.');
      len = (int) strlen (dot+1);
      if (len < 9) tmp = tmp * pow (10.0, 9-len);
      ll = (long) (tmp / 1000000);
      tmp = tmp - ll * 1000000;
      lu = (long) (tmp / 1000);
      la = (long) (tmp - lu * 1000);
    }
    if (!YearWithin(ly) && !ValidateYMD(ly,lm,ld)) return ILLEGAL_TT2000_VALUE;
    return CDF_TT2000_from_UTC_parts ((double) ly, (double) lm, (double) ld,
                                      (double) lh, (double) ln, (double) ls,
                                      (double) ll, (double) lu, (double) la);
  } else
    return ILLEGAL_TT2000_VALUE;
}

/******************************************************************************
* AddOpenCDFCount.
******************************************************************************/

VISIBLE_PREFIX void AddOpenCDFsCount ()
{
  ++openCDF64s;
}

/******************************************************************************
* ReduceOpenCDFCount.
******************************************************************************/

VISIBLE_PREFIX void ReduceOpenCDFsCount ()
{
  --openCDF64s;
  if (openCDF64s == 0) CDFClearLeapSecondsTable();
}

