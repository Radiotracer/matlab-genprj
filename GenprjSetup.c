/*
	GenprjSetup.c

	$Id: GenprjSetup.c 9 2006-09-19 20:17:12Z mjs $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <mip/miputil.h>
#include <mip/printmsg.h>
#include <mip/errdefs.h>
#include <mip/getparms.h>
#include <mip/getline.h>
#include <mip/irl.h>

#include "protos.h"

char *pchGetNormBase(char *pchBase)
{
	char *pchTmpDir;
	char *pchNormBase;
	int bNormInMemory;
	int bFound;

	pchTmpDir=pchIrlStrdup(pchGetStrParm("tmpdir",&bFound,"/var/tmp"));
	bNormInMemory = bGetBoolParm("norm_in_memory", &bFound, FALSE);
	if (bNormInMemory)
		return NULL;
	pchNormBase = pvIrlMalloc(strlen(pchTmpDir) + strlen(pchBase)+2,
		"GetNormBase:pchNormBase");
	sprintf(pchNormBase,"%s/%s",pchTmpDir,pchBase);
	return pchNormBase;
}

void vGetPrjParms(IrlParms_t *psParms, Options_t *psOptions)
/* get Parameters that are needed for reconstruction*/
{
	int bFound, bFoundAll, bDrfFromFile;
	int i;

	vPrintMsg(4,"GetPrjParms\n");
	psParms->BinWidth = (float)dGetDblParm("pixwidth",&bFound, 0.0);
   psParms->NumViews = iGetIntParm("nang", &bFound, 1);
	psParms->fScatEstFac=dGetDblParm("scat_est_fac", &bFound, 1.0);
	i=psParms->SrfCollapseFac=iGetIntParm("srf_collapse_fac", &bFound, 1);
	if (i != 1 && i != 2 && i != 4){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE,
							"SrfSetup",
							"SrfFactor must be 1 (no collapse), 2 or 4",
							i);
	}

	psParms->fAtnScaleFac=dGetDblParm("AtnMapFac", &bFound, 1.0);
	bDrfFromFile = bGetBoolParm("drf_from_file",&bFound, FALSE);
	if (psOptions->bModelDrf && (!bDrfFromFile)){
		psParms->fHoleLen=dGetDblParm("collthickness", &bFound, 1.0);
		bFoundAll=bFound;
		psParms->fHoleDiam=dGetDblParm("holediam", &bFound, 1.0);
		bFoundAll &= bFound;
		psParms->fBackToDet=dGetDblParm("gap", &bFound, 1.0);
		bFoundAll &= bFound;
		psParms->fIntrinsicFWHM=dGetDblParm("intrinsicFWHM", &bFound, 1.0);
		bFoundAll &= bFound;
		if (!bFoundAll)
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vGetPrjParms",
	"Missing gap, collthickness holediam or intrinsicfwhm with DRF modeling");
	}else{
		psParms->fHoleLen=0.0;
		psParms->fHoleDiam=0.0;
		psParms->fBackToDet=0.0;
		psParms->fIntrinsicFWHM=0.0;
	}
	if (psOptions->bModelDrf){
		psOptions->bFFTConvolve=bGetBoolParm("fft_convolve", &bFound, FALSE);
	}else{
		psOptions->bFFTConvolve=FALSE;
	}
	psOptions->fAtnMapThresh=dGetDblParm("atnmap_support_thresh",&bFound, 0.0);
	psOptions->bUseContourSupport = bGetBoolParm("use_contour_support",&bFound,0);
	psOptions->fAtnMapThresh *= psParms->BinWidth;
	psOptions->iMsgLevel=iGetIntParm("debug_level",&bFound, 4);
	psOptions->iAxialPadLength=iGetIntParm("axial_pad_length",&bFound,0);
	psOptions->iAxialAvgLength=iGetIntParm("axial_avg_length",&bFound,0);
}

int iParseModelString(char *pchModelStr, char *pchName)
/* returns the integer effects string corresponding to pchModelStr.
	pchModel string must contain only characters a, d, or s. Case is
	not significant. If a is found then MODEL_ATN is set. If d is found
	then MODEL_DRF is set. If "s" is found then MODEL_SRF is set. These
	MODEL_STR bits are all or-ed together and returned. The pchName is the
	name of the model being parsed and is used in reporting errors
*/
{
	int iModel=0;

	while ((*pchModelStr) != '\0'){
		switch (tolower(*pchModelStr)){
			case 'a': 
				iModel |= MODEL_ATN;
				break;
			case 'd':
				iModel |= MODEL_DRF;
				break;
			case 's':
				iModel |= MODEL_SRF;
				break;
			default:
				vErrorHandler(ECLASS_WARN, ETYPE_ILLEGAL_VALUE, "ParseModelString",
					"illegal character (%c) found in %s. Ignoring.",
					*pchModelStr, pchName);
				break;
		}
		pchModelStr++;
	}
	return iModel;
}

void vGetEffectsToModel(int *pbModelAtn, int *pbModelDrf, int *pbModelSrf)
/* Get Effects to model (e.g. atten, drf, srf)*/
{
	char *pch;
	int bFound;
	int iModel;
	
	vPrintMsg(4,"GetEffects\n");
	pch = pchGetStrParm("model",&bFound,"");
	pch = pchGetStrParm("prjmodel",&bFound,pch);
	if (bFound) 
		iModel= iParseModelString(pch, "prjmodel");
	else
		iModel= iParseModelString(pch, "model");
	*pbModelAtn = iModel & MODEL_ATN;
	*pbModelDrf = iModel & MODEL_DRF;
	*pbModelSrf = iModel & MODEL_SRF;
}

PrjView_t *psSetupPrjViews(IrlParms_t *psParms)
/* sets up orbit (defined by radius of rotation and center of rotation in views
	structure. Handles either circular orbit (Single cor and ror) or
	arbitrary orbits (dift cor and ror for each view read from a file) */
{
	float fCFCR;
	float fPwid;
	float fModCntr;
	int bFound, iLine;
	char *pchLine;
	int bNonCircOrbit;
	char *pchOrbitFileName=NULL;
	FILE *fpOrbitFile;
	int iAng;
	PrjView_t *psViews ;
	double dDelAng, dAngleStart, dAngleRange;
	int iNumAngles;

	vPrintMsg(4,"SetupPrjViews\n");

	iNumAngles = psParms->NumViews;
	psViews = pvIrlMalloc(sizeof(PrjView_t)*iNumAngles,"SetupViews:Views");

	dAngleStart = (float)dGetDblParm("angle_start",&bFound, 0);
	dAngleRange = (float)dGetDblParm("angle_range",&bFound, 360);
	dDelAng = dAngleRange / iNumAngles;
	for(iAng = 0; iAng < iNumAngles; ++iAng){
		psViews[iAng].Angle = iAng;
		psViews[iAng].Angle = (dAngleStart + dDelAng*iAng)*M_PI/180.0;
	}

	bNonCircOrbit = bGetBoolParm("Noncircular_orbit",&bFound,FALSE);

	if (! bNonCircOrbit){
  	/* for circular orbit the radius of rotation and center of rotation
	 	 are the same */
		fCFCR = dGetDblParm("Cor2Col", &bFound, 0);
		if (fCFCR < 0.0 )
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, 
				"vSetupOrbit", "Cor2Col must be >= 0.0");
		for (iAng=0; iAng<iNumAngles; iAng++) {
	   	psViews[iAng].CFCR = fCFCR;
		}
		vPrintMsg(7,"Circular Orbit: ror=%.1f cm\n", fCFCR);
	}else{
		vPrintMsg(7,"NonCirular Orbit:\n");
		vPrintMsg(8,"  Angle      Ror      Cor\n");
		pchOrbitFileName = pchGetStrParm("orbit_file",0,NULL);
		fpOrbitFile = fopen(pchOrbitFileName,"rt");
		if(fpOrbitFile == NULL)
	  		 vErrorHandler(ECLASS_FATAL, ETYPE_IO, "vSetupOrbit",
				"cannot open file %s for read",pchOrbitFileName);
		for (iAng=0; iAng<iNumAngles; iAng++){
			pchLine=pchGetLine(fpOrbitFile,0,&iLine);
		  	if (sscanf(pchLine, " %f", &fCFCR) != 1)
				vErrorHandler(ECLASS_FATAL, ETYPE_IO, "vSetupOrbit",
					"Error reading from orbit file %s for angle %d", 
					pchOrbitFileName,iAng);
		   psViews[iAng].CFCR = fCFCR;
			vPrintMsg(8,"%7d %8.1f\n",iAng,psViews[iAng].CFCR);
		}
		fclose(fpOrbitFile);
	}
	
	/* with the assumptions that the pixel and bin widths are the same,
		the number of pixels and number of bins are the same,
		and all the transaxial bins were measured, then Left and Right are
		easy to calculate
	*/
	for (iAng=0; iAng<iNumAngles; iAng++){
		psViews[iAng].Left=0;
		psViews[iAng].Right=psParms->NumPixels-1;
	}
	return psViews;
}

char *pchGetSrfKrnlFname(void)
{
	char *pch;
	int bFound;

	pch=pchGetStrParm("srf_krnl_file",&bFound,"");
	if (!bFound)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetSrfKrnlFname",
		"Scatter modeling requested but srf_krnl_file not in parameter file");
	return pchIrlStrdup(pch);
}

char *pchGetDrfTabFname(void)
{
	char *pch;
	int bFound;

	pch=pchGetStrParm("drf_tab_file",&bFound,"");
	if (!bFound)
		return NULL;
	return pchIrlStrdup(pch);
}


