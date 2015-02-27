/*
	genprj.c
	
	$Id: genprj.c 12 2009-03-06 21:03:13Z mjs $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mip/errdefs.h>
#include <mip/printmsg.h>
#include <mip/cmdline.h>
#include <mip/getparms.h>
#include <mip/miputil.h>
#include <mip/imgio.h>
#include <mip/irl.h>

#include "protos.h"
#include "buildinfo.h"
#include "mex.h"


#define IMAGE_EXTENSION ".im"
#define MAXNAMES 3

struct {
	int iNumBins;
	int iNumSlices;
	char *pchOutBase;
   IMAGE *pOutPrjImage;
   PrjView_t *ppView;
   float *pfUsrPrjData;
   IrlParms_t *psParm;
} sPrjCallbackData;

void *pvLocalMalloc(int iSize, char *pchName)
{
	void *pvPtr;
	pvPtr=malloc(iSize);
	if (pvPtr == NULL){
		fprintf(stderr,"LocalMalloc: unable to allocate %d bytes for %s",
			iSize, pchName);
		exit(1);
	}
	return(pvPtr);
}

int exists(const char *fname)
{
    FILE *file;
    if (file = fopen(fname, "r"))
    {
        fclose(file);
        return 1;
    }
    return 0;
}
		
void vGetmxImageSizesForPrj(mxArray *actImg, IrlParms_t *psParms)
{
	int nDim = mxGetNumberOfDimensions(actImg);
	const mwSize *dims = mxGetDimensions(actImg);

	if (nDim<2 || nDim>3)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vGetmxImageSizesForPrj",
		"Image Dimension = %d. The input activity image should be 2D or 3D", nDim);


	if (dims[0] != dims[1])
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vGetmxImageSizesForPrj",
		"X and Y dimension of activity image should be equal : x-dimension = %d, Y-dimension = %d",
		dims[0], dims[1]);
	psParms->NumPixels = dims[0];

	// The slice selection in the parameter file is disabled.Use all slices.
	if (nDim == 2)  
		psParms->NumSlices = 1;
	else
		psParms->NumSlices = dims[2];

}



float* ToFloatArray(const mxArray* input)
{
	float* output;
	int i;
	mwSize ndim = mxGetNumberOfDimensions(input);
	mwSize *tdims = mxGetDimensions(input);
	unsigned int numel = mxGetNumberOfElements(input);

	output = pvLocalMalloc(numel*sizeof(float), "ToFloatArray:output");
	
	mwSize dims[3];
	dims[0] = tdims[0];
	dims[1] = tdims[1];
	if (ndim == 2)
	{
		dims[2] = 1;
	}
	else if (ndim == 3)
	{
		dims[2] = tdims[2];
	}
	else
		mexErrMsgTxt("Input image matrix must be 2d or 3d.");


	switch (mxGetClassID(input))
	{
	case mxDOUBLE_CLASS:
	{
						   double* inptr = (double*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2];s++)						   
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							  output[s*dims[0] * dims[1] + m*dims[1]+n] = inptr[i];
							  i++;
						   }
						   break;
	}
	case mxSINGLE_CLASS:
	{
						   float *inptr = (float*)mxGetData(input);
						   	i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2];s++)						   
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							  output[s*dims[0] * dims[1] + m*dims[1]+n] = inptr[i];
							  i++;
						   }
						   break;
	}
	case mxINT8_CLASS:
	{
						 signed char* inptr = (signed char*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						 break;
	}
	case mxUINT8_CLASS:
	{
						  unsigned char* inptr = (unsigned char*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						  break;
	}
	case mxINT16_CLASS:
	{
						  short int* inptr = (short int*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						  break;
	}
	case mxUINT16_CLASS:
	{
						   unsigned short int* inptr = (unsigned short int*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						   break;
	}
	case mxINT32_CLASS:
	{
						  int* inptr = (int*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						  break;
	}
	case mxUINT32_CLASS:
	{
						   unsigned int* inptr = (unsigned int*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						   break;
	}

	case mxINT64_CLASS:
	{
						  int64_T* inptr = (int64_T*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						  break;
	}
	case mxUINT64_CLASS:
	{
						   uint64_T* inptr = (uint64_T*)mxGetData(input);
						   i = 0;
						   //Reorder the matrix into row major.
						   for (int s = 0; s < dims[2]; s++)
						   for (int n = 0; n < dims[1]; n++)
						   for (int m = 0; m < dims[0]; m++)
						   {
							   output[s*dims[0] * dims[1] + m*dims[1] + n] = inptr[i];
							   i++;
						   }
						   break;
	}

	default:
		mexErrMsgTxt("The class of input array is not numerical.");
		break;
	}

	return output;
}


// The matlab interface function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
	if (nrhs<2 || nrhs>3 || !mxIsChar(prhs[0]) ||!mxIsNumeric(prhs[1])){
		mexErrMsgTxt("Invalid input. \nUsage:prjimg=genprj('genprj.par',actimg,[atnmap])");
	}		

	IrlParms_t sIrlParms;
	Options_t sOptions;
	PrjView_t *psViews;
	float *pfPrjImage = NULL, *pfAtnMap = NULL,*pfActImage = NULL,*pfScatterEstimate = NULL;
	char *pchDrfTabFile=NULL, *pchSrfKrnlFile=NULL,*pchLogFile=NULL, *pchMsgFile=NULL, *paraFileName = NULL;
	float fPrimaryFac = 1;
	int iMsgLevel = 4;
	int bFound;
	int bModelAtn = 0, bModelSrf = 0, bModelDrf = 0;
	
	long m = mxGetM(prhs[0]); long n = mxGetN(prhs[0]);
	paraFileName = (char*)mxCalloc(n + 1, sizeof(char));
	int status = mxGetString(prhs[0], paraFileName, n + 1);
	if (status != 0)
		mexErrMsgTxt("Failed to get the parameter file name! \n");
	if (!exists(paraFileName))  
		mexErrMsgTxt("Parameter file does not exist! \n");	
	
	vReadParmsFile(paraFileName);
	iMsgLevel = iGetIntParm("debug_level", &bFound, iMsgLevel);
	vSetMsgLevel(iMsgLevel);

	pfActImage = ToFloatArray(prhs[1]);

	vGetEffectsToModel(&bModelAtn, &bModelDrf, &bModelSrf);
	sOptions.bModelDrf = bModelDrf;

	/* Get Sizes of projection and reconstructed images and number of bins,
	slices, and angles. The number of bins, angles and slices are stored
	in the Parms structure.*/
	vGetmxImageSizesForPrj(prhs[1], &sIrlParms);
	vGetPrjParms(&sIrlParms, &sOptions);

	float fTrueBinWidth = (float)dGetDblParm("binwidth", &bFound, sIrlParms.BinWidth);
	if (fTrueBinWidth != sIrlParms.BinWidth){
		vErrorHandler(ECLASS_FATAL, ETYPE_USAGE, "mexFunction", "%s\n",
			"Bin size must equals pixel size!");
	}

	/* Get Attenuation image. Parameters in psIrlParms are used to determine
	how many are needed */
	if (bModelAtn || bModelSrf){
		if (nrhs!=3)
			mexErrMsgTxt("\n The attenuation map must be provided to model attenuation or scatter. \n");
		//Todo: check the size of attenuation map
		pfAtnMap = ToFloatArray(prhs[2]);
	}

	psViews = psSetupPrjViews(&sIrlParms);

	pfPrjImage = pvLocalMalloc(sIrlParms.NumPixels*sIrlParms.NumSlices*sIrlParms.NumViews*sizeof(float),"mexfunction:pfPrjImage");

	if (bModelSrf)
		pchSrfKrnlFile = pchGetSrfKrnlFname();
	else
		pchSrfKrnlFile = NULL;

	/*Scatter Estimate to be added to computed projection data
	This can be done in addition to the model scatter to allow, e.g.,
	for compensation for downscatter*/
	pfScatterEstimate = pfGetScatterEstimate(&sIrlParms, psViews);

	if (bModelDrf)
		pchDrfTabFile = pchGetDrfTabFname();
	else
		pchDrfTabFile = NULL;

	char* pch = pchGetStrParm("msg_file", &bFound, "");
	if (pch == NULL || *pch == '\0')
		pchMsgFile = NULL;
	else
		pchMsgFile = pchIrlStrdup(pch);

	pch = pchGetStrParm("log_file", &bFound, "");
	if (pch == NULL || *pch == '\0')
		pchLogFile = NULL;
	else
		pchLogFile = pchIrlStrdup(pch);

	PrintTimes("Done Initialization");
	iDoneWithParms();	
	
	int err_num = IrlGenprj(&sIrlParms, &sOptions, psViews,
		pchDrfTabFile, pchSrfKrnlFile,
		NULL,pfScatterEstimate, pfActImage, pfAtnMap, pfPrjImage, fPrimaryFac,
		pchLogFile, pchMsgFile);

	if (err_num){
		fprintf(stderr, "fatal error in IrlGenprj: ErrNum=%d\n      %s",
			err_num, pchIrlErrorString());
	}

	// Generate the output projection image
	nlhs = 1;
	mwSize prjdims[3];
	prjdims[0] = sIrlParms.NumSlices;
	prjdims[1] = sIrlParms.NumPixels;
	prjdims[2] = sIrlParms.NumViews; 
	
 	plhs[0] = mxCreateNumericArray(3, prjdims, mxDOUBLE_CLASS, mxREAL);
	double *pOut = mxGetData(plhs[0]);
	int i = 0;
	for (int s = 0; s < prjdims[2]; s++)
	for (int m = 0; m < prjdims[0]; m++)
	for (int n = 0; n < prjdims[1]; n++ )
	{
		pOut[ s*prjdims[0] *prjdims[1]  + n*prjdims[0]+m] = pfPrjImage[i];
		i++;
	}
	 
}

