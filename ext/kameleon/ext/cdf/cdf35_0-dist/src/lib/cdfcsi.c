/******************************************************************************
* Copyright 1996-2013 United States Government as represented by the
* Administrator of the National Aeronautics and Space Administration.
* All Rights Reserved.
******************************************************************************/
/******************************************************************************
*
*  NSSDC/CDF			CDF C interface for non-macro'ed functions
*				of the Standard Interface.
*
*  Version 2.4a, 8-Mar-97, Hughes STX.
*
*  Modification history:
*
*   V1.0   1-Jun-91, J Love	Original version (for CDF V2.1).  This is a
*				combination of cdf.c, cdfattr.c and cdfvar.c.
*				Most of these functions can be replaced by
*				the macros in 'cdf.h'.
*   V1.1  30-Jul-91, J Love	Use 'CDFlib'.
*   V2.0  10-Feb-92, J Love	IBM PC port.
*   V2.1  21-Aug-92, J Love	CDF V2.3 (shareable/NeXT/zVar).
*   V2.2  18-Oct-93, J Love	CDF V2.4.
*   V2.3   9-Nov-94, J Love	CDF V2.5.
*   V2.4  14-Feb-96, J Love	CDF V2.6 (renamed - previously `cdf_c_if.c').
*   V2.4a  8-Mar-97, J Love	Windows NT for Visual C++ 4 on an IBM PC.
*   V3.0  28-Aug-01, M Liu      Add CDFgetrVarsRecordData, CDFgetzVarsRecordData,
*                               CDFputrVarsRecordData, CDFputzVarsRecordData.
*   V3.1  31-May-05, M Liu      Replaced CDFgetrVarsRecordData and
*                               CDFgetzVarsRecordData with a general function.
*                               So, are the CDFputrVarsRecordData and
*                               CDFputzVarsRecordData.
*   V3.2  13-Apr-07, D. Han     Modified CDFattrInquire not to select CDF id
*                               twice (performance improvement).
*   V3.3  22-Oct-08, M. Liu     Modified CDFgetVarsRecordDatabyNames and
*                               CDFputVarsRecordDatabyNames allocate buffers 
*                               based on the number of variables involved.
*   V3.4  03-Jan-12, M. Liu     Added CDFgetVarAllRecordsByVarID,
*                               CDFgetVarRangeRecordsByVarID,
*                               CDFgetVarAllRecordsByVarName,
*                               CDFgetVarRangeRecordsByVarName functions, and
*                               a set of similar functions for put operations.
*   V3.5  28-May-12, M. Liu     Added CDFinsertVarRecordsByVarID and
*                               CDFinsertVarRecordsByVarName.
*
******************************************************************************/

#include "cdflib.h"
#include "cdflib64.h"

static int backward = 0;
static int setFileBackwardFlag = 0;

static long checksum = 0;
static int setChecksumFlag = 0;

/******************************************************************************
* CDFcreateCDF.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFcreateCDF (cdfName, id)
char    *cdfName;	/* In -- CDF name. */
CDFid   *id;            /* Out -- CDF id. */
{
  CDFstatus pStatus = CDF_OK;
  long dimSizes[1] = {0};
  if (!sX(CDFlib(CREATE_, CDF_, cdfName, 0L, dimSizes, id,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFattrInquire.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use. 
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFattrInquire (id, attrNum, attrName, attrScope,
                                         maxEntry)
CDFid   id;             /* In -- CDF id. */
long    attrNum;        /* In -- Attribute number. */
char    *attrName;      /* Out -- Attribute name. */
long    *attrScope;     /* Out -- Attribute scope. */
long    *maxEntry;      /* Out -- Maximum gEntry/rEntry number used. */
{
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, attrScope,
                 NULL_), &pStatus)) return pStatus;
  if (!sX(CDFlib(GET_, ATTR_NAME_, attrName,
                       BOO(GLOBALscope(*attrScope),ATTR_MAXgENTRY_,
                                                   ATTR_MAXrENTRY_), maxEntry,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFinquireAttr.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFinquireAttr (id, attrNum, attrName, scope, 
					 maxgEntry, maxrEntry, maxzEntry)
CDFid   id;             /* In -- CDF id. */
long    attrNum;        /* In -- Attribute number. */
char    *attrName;      /* Out -- Attribute name. */
long    *scope;         /* Out -- Attribute scope. */
long    *maxgEntry;     /* Out -- Maximum gEntry number used for gAttribute. */
long    *maxrEntry;     /* Out -- Maximum rEntry number used for vAttribute. */
long    *maxzEntry;     /* Out -- Maximum zEntry number used for vAttribute. */
{
 
  CDFstatus pStatus = CDF_OK;

  *maxgEntry = -1L;
  *maxrEntry = -1L;
  *maxzEntry = -1L;

  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, scope,
                 NULL_), &pStatus)) return pStatus;
  if (!sX(CDFlib(GET_, ATTR_NAME_, attrName,
                       (BOO(GLOBALscope(*scope),ATTR_MAXgENTRY_,
                                		ATTR_MAXrENTRY_)), 
                       (BOO(GLOBALscope(*scope),maxgEntry,maxrEntry)),
                 NULL_), &pStatus)) return pStatus;
  if (!GLOBALscope(*scope)) {
    if (!sX(CDFlib(GET_, ATTR_MAXzENTRY_, maxzEntry,
		   NULL_), &pStatus)) return pStatus;
  } 

  return pStatus;
}

/******************************************************************************
* CDFinquireAttrEntry.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/
  
VISIBLE_PREFIX CDFstatus CDFinquireAttrEntry (id, grzEntry, attrNum, 
                                              entryNum, dataType, numElems)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    entryNum;       /* In -- gEntry/rEntry/zEntry number. */
long    *dataType;      /* Out -- gEntry/rEntry/zEntry data type. */
long    *numElems;      /* Out -- gEntry/rEntry/zEntry number of elements. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 0) && (grzEntry != 1)) 
    return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry == 1)) return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3),zENTRY_,
                              (BOO(GLOBALscope(scope),gENTRY_,
                                   rENTRY_))), entryNum,
                 GET_, BOO((grzEntry == 3),zENTRY_DATATYPE_,
                           (BOO(GLOBALscope(scope),gENTRY_DATATYPE_,
                                rENTRY_DATATYPE_))), dataType,
                       BOO((grzEntry == 3),zENTRY_NUMELEMS_,
                           (BOO(GLOBALscope(scope),gENTRY_NUMELEMS_,
                                rENTRY_NUMELEMS_))), numElems,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFputAttrEntry.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFputAttrEntry (id, grzEntry, attrNum, entryNum, 
                                          dataType, numElems, value)
CDFid	id;		/* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long	attrNum;	/* In -- Attribute number. */
long	entryNum;	/* In -- gEntry/rEntry/zEntry number. */
long	dataType;	/* In -- gEntry/rEntry/zEntry data type. */
long	numElems;	/* In -- gEntry/rEntry/zEntry number of elements. */
void	*value;		/* In -- Value. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
			  ATTR_, attrNum,
		 GET_, ATTR_SCOPE_, &scope,
		 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 0) && (grzEntry != 1)) 
    return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry == 1)) return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3),zENTRY_,
                              (BOO(GLOBALscope(scope),gENTRY_,
                               rENTRY_))), entryNum,
		 PUT_, BOO((grzEntry == 3),zENTRY_DATA_,
                           (BOO(GLOBALscope(scope),gENTRY_DATA_,
                                rENTRY_DATA_))), dataType, numElems, value,
		 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetAttrEntry.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetAttrEntry (id, grzEntry, attrNum, entryNum,
                                          value)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    entryNum;       /* In -- gEntry/rEntry/zEntry number. */
void    *value;         /* Out -- Value. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 0) && (grzEntry != 1))
    return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry == 1)) return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3),zENTRY_,
                              (BOO(GLOBALscope(scope),gENTRY_,
                               rENTRY_))), entryNum,
                 GET_, BOO((grzEntry == 3),zENTRY_DATA_,
                           (BOO(GLOBALscope(scope),gENTRY_DATA_,
                                rENTRY_DATA_))), value,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFdeleteAttrEntry.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFdeleteAttrEntry (id, grzEntry, attrNum, entryNum)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    entryNum;       /* In -- gEntry/rEntry/zEntry number. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 1)) return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry != 2) && (grzEntry != 3))
    return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3),zENTRY_,
                              (BOO(GLOBALscope(scope),gENTRY_,
                               rENTRY_))), entryNum,
                 DELETE_, BOO((grzEntry == 3),zENTRY_,
                              (BOO(GLOBALscope(scope),gENTRY_, rENTRY_))), 
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFsetAttrEntryDataSpec.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFsetAttrEntryDataSpec (id, grzEntry, attrNum, 
                                                  entryNum, dataType, numElems)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    entryNum;       /* In -- gEntry/rEntry/zEntry number. */
long    dataType;       /* In -- Data type. */
long    numElems;       /* In -- Number of elements. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 1)) return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry != 2) && (grzEntry != 3))
    return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3),zENTRY_,
                              (BOO(GLOBALscope(scope),gENTRY_,
                               rENTRY_))), entryNum,
                 PUT_, BOO((grzEntry == 3),zENTRY_DATASPEC_,
                           (BOO(GLOBALscope(scope),gENTRY_DATASPEC_,
                                                   rENTRY_DATASPEC_))), 
                       dataType, numElems,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetAttrNum.
* Can't implement with macro since it is the attribute number which is to be
* returned (unless an error).
******************************************************************************/

VISIBLE_PREFIX long CDFgetAttrNum (id,attrName)
CDFid	id;		/* In -- CDF id. */
char	*attrName;	/* In -- attribute name. */
{
  CDFstatus status;
  long attrNum;
  status = CDFlib (SELECT_, CDF_, id,
		   GET_, ATTR_NUMBER_, attrName, &attrNum,
		   NULL_);
  if (StatusOK(status))
    return attrNum;
  else
    return status;
}

/******************************************************************************
* CDFgetVarNum.
* Can't implement with macro since it is the variable number which is to be
* returned (unless an error).
******************************************************************************/

VISIBLE_PREFIX long CDFgetVarNum (id,varName)
CDFid	id;		/* In -- CDF id. */
char	*varName;	/* In -- variable name. */
{
  CDFstatus status;
  long varNum;
  status = CDFlib (SELECT_, CDF_, id,
		   GET_, zVAR_NUMBER_, varName, &varNum,
		   NULL_);
  if (StatusOK(status))
    return varNum;
  else {
    status = CDFlib (SELECT_, CDF_, id,
	             GET_, rVAR_NUMBER_, varName, &varNum,
		     NULL_);
    if (StatusOK(status))
      return varNum;
  }  
    return status;
}

/******************************************************************************
* CDFgetNumAttrEntries.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use. 
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetNumAttrEntries (id, grzEntry, attrNum, 
                                               numEntries)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    *numEntries;    /* Out -- number of gEntries/rEntries/zEntries. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 1)) return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry != 2) && (grzEntry != 3))
    return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(GET_, BOO((grzEntry == 3), ATTR_NUMzENTRIES_, 
                           (BOO(GLOBALscope(scope), ATTR_NUMgENTRIES_,
                                ATTR_NUMrENTRIES_))), numEntries,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetAttrMaxEntry.  
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/
  
VISIBLE_PREFIX CDFstatus CDFgetAttrMaxEntry (id, grzEntry, attrNum, maxEntry)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    *maxEntry;      /* Out -- Max. number of gEntry/rEntry/zEntry. */
{                
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 1)) return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry != 2) && (grzEntry != 3))
    return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(GET_, BOO((grzEntry == 3), ATTR_MAXzENTRY_,
                           (BOO(GLOBALscope(scope), ATTR_MAXgENTRY_, 
                            ATTR_MAXrENTRY_))), maxEntry,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetAttrEntryDataType.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetAttrEntryDataType (id, grzEntry, attrNum, 
                                                  entryNum, dataType)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* In -- Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    entryNum;       /* In -- gEntry/rEntry/zEntry number. */
long    *dataType;      /* Out -- gEntry/rEntry/zEntry dataType. */
{
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 1)) return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry != 2) && (grzEntry != 3))
    return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3), zENTRY_, 
                              BOO(GLOBALscope(scope), gENTRY_,
                                  rENTRY_)), entryNum,
                 GET_, BOO((grzEntry == 3), zENTRY_DATATYPE_,
                           BOO(GLOBALscope(scope), gENTRY_DATATYPE_,
                               rENTRY_DATATYPE_)), dataType,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetAttrEntryNumElements.
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetAttrEntryNumElements (id, grzEntry, attrNum, 
                                                     entryNum, numElements)
CDFid   id;             /* In -- CDF id. */
int     grzEntry;       /* Flag for g/r/zEntry. */
long    attrNum;        /* In -- Attribute number. */
long    entryNum;       /* In -- gEntry/rEntry/zEntry number. */
long    *numElements;   /* Out -- gEntry/rEntry/zEntry numElements. */
{                
  long scope;
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,         
                 GET_, ATTR_SCOPE_, &scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(scope) && (grzEntry != 1)) return ILLEGAL_FOR_SCOPE;
  if (!GLOBALscope(scope) && (grzEntry != 2) && (grzEntry != 3))
    return ILLEGAL_FOR_SCOPE;

  if (!sX(CDFlib(SELECT_, BOO((grzEntry == 3), zENTRY_, 
                              BOO(GLOBALscope(scope), gENTRY_,
                                  rENTRY_)), entryNum,
                 GET_, BOO((grzEntry == 3), zENTRY_NUMELEMS_, 
                           BOO(GLOBALscope(scope), gENTRY_NUMELEMS_,
                               rENTRY_NUMELEMS_)), numElements,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;         
}                
                 
/******************************************************************************
* CDFgetVarRecordData.  
* Acquire a full record data from a given record of a variable.
* Retrieved data are filled into the buffer.
******************************************************************************/
  
VISIBLE_PREFIX CDFstatus CDFgetVarRecordData (id, zVar, varNum, recNum, buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
long    recNum;         /* In -- Record number to read. */
void    *buffer;        /* Out -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  long numVars, varNums[CDF_MAX_DIMS];

  numVars = 1;
  varNums[0] = varNum;

  if (!sX(CDFlib(SELECT_, CDF_, id,
                          BOO((zVar == 1),zVARs_RECNUMBER_,rVARs_RECNUMBER_),
                          recNum,
                 GET_, BOO((zVar == 1),zVARs_RECDATA_,rVARs_RECDATA_), numVars,
                       varNums, buffer,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFputVarRecordData.
* Write a full record data for a given record of a variable.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFputVarRecordData (id, zVar, varNum, recNum, buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
long    recNum;         /* In -- Record number to write. */
void    *buffer;        /* In -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  long numVars, varNums[CDF_MAX_DIMS];

  numVars = 1;
  varNums[0] = varNum;

  if (!sX(CDFlib(SELECT_, CDF_, id,
                          BOO((zVar == 1),zVARs_RECNUMBER_,rVARs_RECNUMBER_),
                          recNum,
                 PUT_, BOO((zVar == 1),zVARs_RECDATA_,rVARs_RECDATA_), numVars,
                       varNums, buffer,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetVarAllRecordsByVarID.  
* Acquire data from all records for a given variable (by zVar flag and id).
* Retrieved data are filled into the buffer.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetVarAllRecordsByVarID (id, zVar, varNum, buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
void    *buffer;        /* Out -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  int i;
  long lastRec, numDims;
  long indices[CDF_MAX_DIMS], intervals[CDF_MAX_DIMS], dimSizes[CDF_MAX_DIMS];

  if (zVar == 1) {
    if (!sX(CDFgetzVarNumDims(id,varNum,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetzVarDimSizes(id,varNum,dimSizes), &pStatus)) return pStatus;
    if (!sX(CDFgetVarMaxWrittenRecNum(id,1,varNum,&lastRec), &pStatus))
      return pStatus;
  } else {
    if (!sX(CDFgetrVarsNumDims(id,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetrVarsDimSizes(id,dimSizes), &pStatus)) return pStatus;
    if (!sX(CDFgetVarMaxWrittenRecNum(id,0,varNum,&lastRec), &pStatus))
      return pStatus;
  }
  for (i = 0; i < (int) numDims; ++i) {
    indices[i] = 0L;
    intervals[i] = 1L;
  }
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          BOO((zVar == 1),zVAR_,rVAR_), varNum,
                          BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),0L,
                          BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                          lastRec+1,
                          BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                          1L,
                          BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                          indices,
                          BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                          dimSizes,
                          BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                          intervals,
                 GET_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),buffer,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetVarAllRecordsByVarName.  
* Acquire data from all records for a given variable (by its unique name).
* Retrieved data are filled into the buffer.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetVarAllRecordsByVarName (id, varName, buffer)
CDFid   id;             /* In -- CDF id. */
char    *varName;       /* In -- Variable name. */
void    *buffer;        /* Out -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  int zVar;
  long varNum;

  varNum = CDFgetVarNum(id, varName);
  if (varNum < 0) return NO_SUCH_VAR;
  if (CDFconfirmVarExistence(id, 1, varName) == CDF_OK) zVar = 1;
  else zVar = 0;
  if (!sX(CDFgetVarAllRecordsByVarID(id,zVar,varNum,buffer), &pStatus))
    return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFputVarAllRecordsByVarID.
* Write a numbr of records for a given variable (by zVar flag and id), 
* staring from record 0.
******************************************************************************/
  
VISIBLE_PREFIX CDFstatus CDFputVarAllRecordsByVarID (id, zVar, varNum, numRecs,
                                                     buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
long    numRecs;        /* In -- Number of Records to write. */
void    *buffer;        /* In -- Buffer for holding data. */ 
  
{ 
  CDFstatus pStatus = CDF_OK;
  int i;
  long numDims, maxRec;
  long indices[CDF_MAX_DIMS], intervals[CDF_MAX_DIMS], dimSizes[CDF_MAX_DIMS];
  
  if (zVar == 1) { 
    if (!sX(CDFgetzVarNumDims(id,varNum,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetzVarDimSizes(id,varNum,dimSizes), &pStatus)) return pStatus;
  } else {
    if (!sX(CDFgetrVarsNumDims(id,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetrVarsDimSizes(id,dimSizes), &pStatus)) return pStatus;
  } 
  for (i = 0; i < (int) numDims; ++i) {
    indices[i] = 0L;
    intervals[i] = 1L;
  }                       
  if (!sX(CDFgetVarMaxWrittenRecNum(id,zVar,varNum,&maxRec),&pStatus))
    return pStatus;
  if (!sX(CDFlib(SELECT_, CDF_, id, 
                          BOO((zVar == 1),zVAR_,rVAR_), varNum,
                          BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),0L,
                          BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                          numRecs,
                          BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                          1L,
                          BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                          indices,
                          BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                          dimSizes,
                          BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                          intervals,
                 PUT_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),buffer,
                 NULL_), &pStatus)) return pStatus;
  /****************************************************************************
  * Update maximum record numbers.
  ****************************************************************************/
  if (maxRec > (numRecs-1)) {
    if (!sX(CDFlib(DELETE_, BOO((zVar == 1),zVAR_RECORDS_,rVAR_RECORDS_),
                            numRecs, maxRec,
                   NULL_), &pStatus)) return pStatus;
  }
  return pStatus;
}

/******************************************************************************
* CDFputVarAllRecordsByVarName.
* Write a number of records for a given variable (by its unique name), staring
* from record 0.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFputVarAllRecordsByVarName (id, varName, numRecs,
                                                       buffer)
CDFid   id;             /* In -- CDF id. */
char    *varName;       /* In -- Variable name. */
long    numRecs;        /* In -- Number of Records to write. */
void    *buffer;        /* In -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  int zVar;
  long varNum;

  varNum = CDFgetVarNum(id, varName);
  if (varNum < 0) return NO_SUCH_VAR;
  if (CDFconfirmVarExistence(id, 1, varName) == CDF_OK) zVar = 1;
  else zVar = 0;
  if (!sX(CDFputVarAllRecordsByVarID(id,zVar,varNum,numRecs,buffer), &pStatus))
    return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetVarRangeRecordsByVarID.  
* Acquire data from a range of records for a given variable (by zVar flag and
* id). Retrieved data are filled into the buffer.      
******************************************************************************/
  
VISIBLE_PREFIX CDFstatus CDFgetVarRangeRecordsByVarID (id, zVar, varNum,
                                                       startRec, stopRec,
                                                       buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
long    startRec;       /* In -- Starting record number. */
long    stopRec;        /* In -- Ending record number. */
void    *buffer;        /* Out -- Buffer for holding data. */
  
{ 
  CDFstatus pStatus = CDF_OK;
  int i;
  long numDims;
  long indices[CDF_MAX_DIMS], intervals[CDF_MAX_DIMS], dimSizes[CDF_MAX_DIMS];
    
  if (zVar == 1) {
    if (!sX(CDFgetzVarNumDims(id,varNum,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetzVarDimSizes(id,varNum,dimSizes), &pStatus)) return pStatus;
  } else {
    if (!sX(CDFgetrVarsNumDims(id,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetrVarsDimSizes(id,dimSizes), &pStatus)) return pStatus;
  } 
  for (i = 0; i < (int) numDims; ++i) {
    indices[i] = 0L;
    intervals[i] = 1L;    
  }                       
  if (!sX(CDFlib(SELECT_, CDF_, id, 
                          BOO((zVar == 1),zVAR_,rVAR_), varNum,
                          BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),
                          startRec,
                          BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                          stopRec-startRec+1,
                          BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                          1L,
                          BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                          indices,
                          BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                          dimSizes,
                          BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                          intervals,
                 GET_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),buffer,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetVarRangeRecordsByVarName.  
* Acquire data from a range of records for a given variable (by its unique 
* name). Retrieved data are filled into the buffer.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetVarRangeRecordsByVarName (id, varName,
                                                         startRec, stopRec,
                                                         buffer)
CDFid   id;             /* In -- CDF id. */
char    *varName;       /* In -- Variable name. */
long    startRec;       /* In -- Starting record number. */
long    stopRec;        /* In -- Ending record number. */
void    *buffer;        /* Out -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  int  zVar;
  long varNum;

  varNum = CDFgetVarNum(id, varName);
  if (varNum < 0) return NO_SUCH_VAR;
  if (CDFconfirmVarExistence(id, 1, varName) == CDF_OK) zVar = 1;
  else zVar = 0;
  if (!sX(CDFgetVarRangeRecordsByVarID(id,zVar,varNum,startRec,stopRec,
                                       buffer), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFputVarRangeRecordsByVarID.                          
* Write a range of records for a given variable (by zVar flag and id).
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFputVarRangeRecordsByVarID (id, zVar, varNum,
                                                       startRec, stopRec,
                                                       buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
long    startRec;       /* In -- Starting record number to write. */
long    stopRec;        /* In -- Ending record number to write. */
void    *buffer;        /* In -- Buffer for holding data. */
  
{ 
  CDFstatus pStatus = CDF_OK;
  int i;
  long numDims;
  long indices[CDF_MAX_DIMS], intervals[CDF_MAX_DIMS], dimSizes[CDF_MAX_DIMS];
    
  if (zVar == 1) {
    if (!sX(CDFgetzVarNumDims(id,varNum,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetzVarDimSizes(id,varNum,dimSizes), &pStatus)) return pStatus;
  } else {
    if (!sX(CDFgetrVarsNumDims(id,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetrVarsDimSizes(id,dimSizes), &pStatus)) return pStatus;
  } 
  for (i = 0; i < (int) numDims; ++i) {
    indices[i] = 0L;
    intervals[i] = 1L;    
  }                       
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          BOO((zVar == 1),zVAR_,rVAR_), varNum,
                          BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),
                          startRec, 
                          BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                          stopRec-startRec+1,
                          BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                          1L,
                          BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                          indices,
                          BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                          dimSizes,
                          BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                          intervals,
                 PUT_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),buffer,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFputVarRangeRecordsByVarName.
* Write a range of records for a given variable (by its unique name).
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFputVarRangeRecordsByVarName (id, varName,
                                                         startRec, stopRec,
                                                         buffer)
CDFid   id;             /* In -- CDF id. */
char    *varName;       /* In -- Variable name. */
long    startRec;       /* In -- Starting record number to write. */
long    stopRec;        /* In -- Ending record number to write. */
void    *buffer;        /* In -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  int zVar;
  long varNum;

  varNum = CDFgetVarNum(id, varName);
  if (varNum < 0) return NO_SUCH_VAR;
  if (CDFconfirmVarExistence(id, 1, varName) == CDF_OK) zVar = 1;
  else zVar = 0;
  if (!sX(CDFputVarRangeRecordsByVarID(id,zVar,varNum,startRec,stopRec,buffer),
                                       &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFinsertVarRecordsByVarID.                          
* Insert a number of records for a given variable (by zVar flag and id).
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFinsertVarRecordsByVarID (id, zVar, varNum,
                                                     startRec, numRecs,
                                                     buffer)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    varNum;         /* In -- Variable number. */
long    startRec;       /* In -- Starting record number to insert. */
long    numRecs;        /* In -- The number of records to insert. */
void    *buffer;        /* In -- Buffer for holding data. */
  
{ 
  CDFstatus pStatus = CDF_OK;
  int i;
  long numDims, lastRec, dataType;
  long indices[CDF_MAX_DIMS], intervals[CDF_MAX_DIMS], dimSizes[CDF_MAX_DIMS];
  long sparseness;

  if (numRecs < 1) return CDF_OK;
  if (zVar == 1) {
    if (!sX(CDFgetVarSparseRecords(id,1,varNum,&sparseness), &pStatus))
      return pStatus;
    if (sparseness != NO_SPARSERECORDS) return CANNOT_INSERT_RECORDS;
    if (!sX(CDFgetzVarNumDims(id,varNum,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetVarDataType(id,1,varNum,&dataType), &pStatus)) return pStatus;
    if (!sX(CDFgetzVarDimSizes(id,varNum,dimSizes), &pStatus)) return pStatus;
    if (!sX(CDFgetVarMaxWrittenRecNum(id,1,varNum,&lastRec), &pStatus))
      return pStatus;
  } else {
    if (!sX(CDFgetVarSparseRecords(id,0,varNum,&sparseness), &pStatus))
      return pStatus;
    if (sparseness != NO_SPARSERECORDS) return CANNOT_INSERT_RECORDS;
    if (!sX(CDFgetrVarsNumDims(id,&numDims), &pStatus)) return pStatus;
    if (!sX(CDFgetVarDataType(id,0,varNum,&dataType), &pStatus)) return pStatus;
    if (!sX(CDFgetrVarsDimSizes(id,dimSizes), &pStatus)) return pStatus;
    if (!sX(CDFgetVarMaxWrittenRecNum(id,0,varNum,&lastRec), &pStatus))
      return pStatus;
  }
  for (i = 0; i < (int) numDims; ++i) {
    indices[i] = 0L;
    intervals[i] = 1L;
  }
  if (startRec <= lastRec) {
    void *tmp;
    long nRecordValues, startRecX, startRecY, recs, recs2;
    int dimN;
    for (dimN = 0, nRecordValues = 1; dimN < numDims; dimN++) {
      nRecordValues *= dimSizes[dimN];
    }
    recs2 = recs = lastRec - startRec + 1;
    while (recs > 0) {
      tmp = cdf_AllocateMemory ((size_t)recs * nRecordValues *
                                CDFelemSize(dataType),
                                NULL);
      if (tmp != NULL) break;
      recs = (int) (0.5 * recs + 0.5); 
    }
    if (tmp == NULL) return BAD_MALLOC;
    for (startRecX = lastRec - recs + 1; recs2 > 0;
         recs2 -= recs, startRecX -= recs) {
      if (startRecX < startRec) {
        recs = recs - (startRec - startRecX);
        startRecX = startRec;
      }
      if (!sX(CDFlib(SELECT_, CDF_, id,
                              BOO((zVar == 1),zVAR_,rVAR_), varNum,
                              BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),
                              startRecX, 
                              BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                              recs,
                              BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                              1L,
                              BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                              indices,
                              BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                              dimSizes,
                              BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                              intervals,
                     GET_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),tmp,
                     NULL_), &pStatus)) return pStatus;
      startRecY = startRecX + numRecs;
      if (!sX(CDFlib(SELECT_, CDF_, id,
                              BOO((zVar == 1),zVAR_,rVAR_), varNum,
                              BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),
                              startRecY, 
                              BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                              recs,
                              BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                              1L,
                              BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                              indices,
                              BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                              dimSizes,
                              BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                              intervals,
                     PUT_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),tmp,
                     NULL_), &pStatus)) return pStatus;
      
    }
    cdf_FreeMemory (tmp, NULL);
  }
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          BOO((zVar == 1),zVAR_,rVAR_), varNum,
                          BOO((zVar == 1),zVAR_RECNUMBER_,rVARs_RECNUMBER_),
                          startRec, 
                          BOO((zVar == 1),zVAR_RECCOUNT_,rVARs_RECCOUNT_),
                          numRecs,
                          BOO((zVar == 1),zVAR_RECINTERVAL_,rVARs_RECINTERVAL_),
                          1L,
                          BOO((zVar == 1),zVAR_DIMINDICES_,rVARs_DIMINDICES_),
                          indices,
                          BOO((zVar == 1),zVAR_DIMCOUNTS_,rVARs_DIMCOUNTS_),
                          dimSizes,
                          BOO((zVar == 1),zVAR_DIMINTERVALS_,rVARs_DIMINTERVALS_),
                          intervals,
                 PUT_, BOO((zVar == 1),zVAR_HYPERDATA_,rVAR_HYPERDATA_),buffer,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFinsertVarRecordsByVarName.
* Insert a number of records for a given variable (by its unique name).
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFinsertVarRecordsByVarName (id, varName,
                                                       startRec, numRecs,
                                                       buffer)
CDFid   id;             /* In -- CDF id. */
char    *varName;       /* In -- Variable name. */
long    startRec;       /* In -- Starting record number to insert. */
long    numRecs;        /* In -- The number of records to insert. */
void    *buffer;        /* In -- Buffer for holding data. */

{
  CDFstatus pStatus = CDF_OK;
  int  zVar;
  long varNum;

  varNum = CDFgetVarNum(id, varName);
  if (varNum < 0) return NO_SUCH_VAR;
  if (CDFconfirmVarExistence(id, 1, varName) == CDF_OK) zVar = 1;
  else zVar = 0;
  if (!sX(CDFinsertVarRecordsByVarID(id,zVar,varNum,startRec,numRecs,buffer),
                                     &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFgetVarsRecordDatabyNames. 
* Acquire a full record data for a given record for a set of the selected 
* r/zVariables. Retrieved data are filled into the buffers that are pointed to
* by the passed array of pointers. The selected variables are identified by 
* their names.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFgetVarsRecordDatabyNames (id, zVar, numVars, 
                                                      varNames, recNum, 
                                                      buffptr)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    numVars;        /* In -- Number of variables. */
char    *varNames[];    /* In -- Array of variable names. */
long    recNum;         /* In -- Record number to read. */ 
void    *buffptr[];     /* Out -- Array of buffer pointers for holding data. */

{ 
  CDFstatus pStatus = CDF_OK;
  long dataType, numElems, numDims; 
  long dimSizes[CDF_MAX_DIMS], dimVarys[CDF_MAX_DIMS];
  long dataTypeSize, recNumValues, *phyRecSize;
  long totalSize, offset, *varNums;
  int i, j;
  void *buffer;

  if (numVars <= 0) return pStatus;

  if (!sX(CDFlib(SELECT_, CDF_, id, 
                 NULL_), &pStatus)) return pStatus;

  totalSize = 0;
  phyRecSize = (long *) cdf_AllocateMemory ((size_t)numVars*sizeof(long)*2, NULL);
  if (phyRecSize == NULL) return BAD_MALLOC;
  varNums = phyRecSize+numVars;

  for (i = 0; i < numVars; i++) {
    if (!sX(CDFlib(GET_, BOO((zVar == 1),zVAR_NUMBER_,rVAR_NUMBER_), 
                         varNames[i], &varNums[i], 
                   NULL_), &pStatus)) return pStatus;
    if (!sX(CDFlib(SELECT_, BOO((zVar == 1),zVAR_,rVAR_), varNums[i],
                   GET_, BOO((zVar == 1),zVAR_DATATYPE_,rVAR_DATATYPE_), 
                         &dataType,
                         BOO((zVar == 1),zVAR_NUMELEMS_,rVAR_NUMELEMS_), 
                         &numElems,
                         BOO((zVar == 1),zVAR_NUMDIMS_,rVARs_NUMDIMS_), 
                         &numDims,
                         BOO((zVar == 1),zVAR_DIMSIZES_,rVARs_DIMSIZES_), 
                         dimSizes,
                         BOO((zVar == 1),zVAR_DIMVARYS_,rVAR_DIMVARYS_), 
                         dimVarys,
                   NULL_), &pStatus)) return pStatus;
    if (!sX(CDFlib(GET_, DATATYPE_SIZE_, dataType, &dataTypeSize,
                   NULL_), &pStatus)) return pStatus;

    if (numDims == 0) {
      numDims = 1;
      dimSizes[0] = 1;
      dimVarys[0] = 0;
    }

    recNumValues = 1;
    for (j = 0; j < numDims; j++)
      if (dimVarys[j]) recNumValues *= dimSizes[j];
    phyRecSize[i] = recNumValues * dataTypeSize * numElems;
    totalSize += phyRecSize[i];
  }

  buffer = (void *) cdf_AllocateMemory ((size_t)totalSize, NULL);
  if (buffer == NULL) return BAD_MALLOC;

  if (!sX(CDFlib(SELECT_, BOO((zVar == 1),zVARs_RECNUMBER_,rVARs_RECNUMBER_), 
                         recNum,
                 GET_, BOO((zVar == 1),zVARs_RECDATA_,rVARs_RECDATA_), numVars,
                       varNums, buffer,
                 NULL_), &pStatus)) return pStatus;

  offset = 0;
  for (i = 0; i < numVars; i++) {
    memcpy((char *) buffptr[i], (char *) buffer+offset, phyRecSize[i]);
    offset += phyRecSize[i];
  }
  cdf_FreeMemory (buffer, NULL);
  cdf_FreeMemory (phyRecSize, NULL);
  return CDF_OK;

}

/******************************************************************************
* CDFputVarsRecordDatabyNames.  
* Write a full record data for a given record for a set of the selected 
* r/zVariables. Input data are pointed to by the passed array of pointers.
******************************************************************************/
                   
VISIBLE_PREFIX CDFstatus CDFputVarsRecordDatabyNames (id, zVar, numVars, 
                                                      varNames, recNum, 
                                                      buffptr)
CDFid   id;             /* In -- CDF id. */
int     zVar;           /* In -- Flag for zVariable. */
long    numVars;        /* In -- Number of variables. */
char    *varNames[];    /* In -- Array of variable names. */
long    recNum;         /* In -- Record number to read. */
void    *buffptr[];     /* In -- Array of buffer pointers for holding data. */
    
{   
  CDFstatus pStatus = CDF_OK;
  long dataType, numElems, numDims; 
  long dimSizes[CDF_MAX_DIMS], dimVarys[CDF_MAX_DIMS];
  long dataTypeSize, recNumValues, *phyRecSize;
  long totalSize, offset, *varNums;
  int i, j;
  void *buffer; 

  if (numVars <= 0) return pStatus;
               
  if (!sX(CDFlib(SELECT_, CDF_, id,
                 NULL_), &pStatus)) return pStatus;

  phyRecSize = (long *) cdf_AllocateMemory ((size_t)numVars*sizeof(long)*2, NULL);
  if (phyRecSize == NULL) return BAD_MALLOC;
  varNums = phyRecSize+numVars;

  totalSize = 0;
  for (i = 0; i < numVars; i++) {
    if (!sX(CDFlib(GET_, BOO((zVar == 1),zVAR_NUMBER_,rVAR_NUMBER_), 
                         varNames[i], &varNums[i],
                   NULL_), &pStatus)) return pStatus;
    if (!sX(CDFlib(SELECT_, BOO((zVar == 1),zVAR_,rVAR_), varNums[i],
                   GET_, BOO((zVar == 1),zVAR_DATATYPE_,rVAR_DATATYPE_), 
                         &dataType,
                         BOO((zVar == 1),zVAR_NUMELEMS_,rVAR_NUMELEMS_), 
                         &numElems,
                         BOO((zVar == 1),zVAR_NUMDIMS_,rVARs_NUMDIMS_), 
                         &numDims,
                         BOO((zVar == 1),zVAR_DIMSIZES_,rVARs_DIMSIZES_), 
                         dimSizes,
                         BOO((zVar == 1),zVAR_DIMVARYS_,rVAR_DIMVARYS_), 
                         dimVarys,
                   NULL_), &pStatus)) return pStatus;
    if (!sX(CDFlib(GET_, DATATYPE_SIZE_, dataType, &dataTypeSize,
                   NULL_), &pStatus)) return pStatus;

    if (numDims == 0) {
      numDims = 1;
      dimSizes[0] = 1;
      dimVarys[0] = 0;
    }

    recNumValues = 1;
    for (j = 0; j < numDims; j++)
      if (dimVarys[j]) recNumValues *= dimSizes[j];
    phyRecSize[i] = recNumValues * dataTypeSize * numElems;
    totalSize += phyRecSize[i];
  }

  buffer = (void *) cdf_AllocateMemory ((size_t)totalSize, NULL);
  if (buffer == NULL) return BAD_MALLOC;

  offset = 0;
  for (i = 0; i < numVars; i++) {
    memcpy((char *) buffer+offset, (char *) buffptr[i], phyRecSize[i]);
    offset += phyRecSize[i];
  }

  if (!sX(CDFlib(SELECT_, BOO((zVar == 1),zVARs_RECNUMBER_,rVARs_RECNUMBER_), 
                         recNum,
                 PUT_, BOO((zVar == 1),zVARs_RECDATA_,rVARs_RECDATA_), numVars,
                       varNums, buffer,
                 NULL_), &pStatus)) return pStatus;

  cdf_FreeMemory (buffer, NULL);
  cdf_FreeMemory (phyRecSize, NULL);
  return CDF_OK;

}

/******************************************************************************
* CDFinquireAttrInfo. (Holding for backward version -- should not use)
* Can't implement with macro because the attribute's scope determines which
* item(s) to use.
******************************************************************************/

VISIBLE_PREFIX CDFstatus CDFinquireAttrInfo (id, zEntry, attrNum, attrName,
                                             scope, maxEntry)
CDFid   id;             /* In -- CDF id. */
int     zEntry;         /* In -- Flag for zEntry. */
long    attrNum;        /* In -- Attribute number. */
char    *attrName;      /* Out -- Attribute name. */
long    *scope;         /* Out -- Attribute scope. */
long    *maxEntry;      /* Out -- Maximum g/r/zEntry number used. */
{
  CDFstatus pStatus = CDF_OK;
  if (!sX(CDFlib(SELECT_, CDF_, id,
                          ATTR_, attrNum,
                 GET_, ATTR_SCOPE_, scope,
                 NULL_), &pStatus)) return pStatus;
  if (GLOBALscope(*scope) && (zEntry == 1)) return ILLEGAL_FOR_SCOPE;
  if (!sX(CDFlib(GET_, ATTR_NAME_, attrName,
                       BOO((zEntry == 1),ATTR_MAXzENTRY_,
                           (BOO(GLOBALscope(*scope),ATTR_MAXgENTRY_,
                                ATTR_MAXrENTRY_))), maxEntry,
                 NULL_), &pStatus)) return pStatus;
  return pStatus;
}

/******************************************************************************
* CDFsetFileBackward. This function has precedence over the environment
* variable in determining whether a backward file is to be created.
******************************************************************************/
  
VISIBLE_PREFIX void CDFsetFileBackward (flag)
int flag;         
{                 
  if (flag != 0) backward = 1;
  else backward = 0;
  /*
   * Set setFileBackward flag. So, if the environment variable is also
   * set, we can ignore that as calling CDFsetFileBackward function has
   * precedence over the environment variable approach.
   */
  setFileBackwardFlag = 1;
} 

/******************************************************************************
* CDFsetFileBackward2. This function is called if the environment variable
* is set to create the backward file. However, if CDFsetFileBackward
* function has been called and as it has precedence over the environment
* variable, calling this function will have no effect at all.
******************************************************************************/

VISIBLE_PREFIX void CDFsetFileBackward2 (flag) 
int flag;
{ 
  if (setFileBackwardFlag == 0) { /* True if CDFsetFileBackward is not */
                                  /* called.                           */
    if (flag != 0) backward = 1;
    else backward = 0;
  }
}

/******************************************************************************
* CDFsetChecksumMode. This function has precedence over the environment
* variable in determining whether a checksum is to be used.
******************************************************************************/

VISIBLE_PREFIX void CDFsetChecksumMode (flag)
long flag;
{
  checksum = flag;
  /*
   * Set setChecksum flag. So, if the environment variable is also
   * set, we can ignore that as calling CDFsetChecksum function has
   * precedence over the environment variable approach.
   */
  setChecksumFlag = 1;
}

/******************************************************************************
* CDFgetFileBackward.
*   Acquires the backward flag defined from CDFsetFileBackward function or
*   the environment variable.
******************************************************************************/

VISIBLE_PREFIX int CDFgetFileBackward ()
{
  return backward;
}

/******************************************************************************
* CDFgetChecksumMode.
*   Acquires the checksum flag defined from CDFsetChecksum function or
*   the environment variable.
******************************************************************************/

VISIBLE_PREFIX long CDFgetChecksumMode ()
{
  return checksum;
}

/******************************************************************************
* CDFsetFileBackwardFlag. (Holding for backward version -- should not use.)
******************************************************************************/

VISIBLE_PREFIX void CDFsetFileBackwardFlag (flag)
int flag;
{
/*
  long version;
  CDFstatus status;
  status = CDFlib (GET_, LIB_VERSION_, &version,
                   NULL_);
  if (version == 3) backward = flag;
*/
  backward = flag;
}

/******************************************************************************
* CDFgetFileBackwardFlag. (Holding for backward version -- should not use.)
******************************************************************************/

VISIBLE_PREFIX int CDFgetFileBackwardFlag ()
{
  if (backward == 0) return 0;
  else return 1;
}

