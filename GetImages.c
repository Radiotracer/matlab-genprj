/**
	@file GetImages.c

	@code $Id: GetImages.c 7 2005-03-23 19:17:08Z mjs duyong$ @endcode

	@brief Several functions to read in the image data required for
	reconstruction.

	@note Define DEBUG while compiling to print out some extra information
	and save images of the modified projections and atten map in reordered
	format (including scaling). See reorder comment below.
	
	If you want these functions to reorder the image pixels to the IRL
	internal format then define REORDER_PIXELS while compiling.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mip/irl.h>
#include <mip/miputil.h>
#include <mip/imgio.h>
#include <mip/errdefs.h>
#include <mip/getparms.h>
#include <mip/printmsg.h>

#include "protos.h"

/**
	 @brief Gets sizes and number of angles for projection data and
	 actitity image.

	 These are obtained from the parameter file and the activity image.
	 
	 @param *pchActImageName - ptr to the activity image data filename.
	 @param *psParms - ptr to the IrlParms_t struct.
	 
	 @return void
*/

void vGetImageSizesForPrj(char *pchActImageName, IrlParms_t *psParms)
{
	int iXdim, iYdim, iZdim;
   int iSliceStart, iSliceEnd;
	int bFound;
	int iNumSlices;
	IMAGE *pActImage;

   int iX, iY, iZ,i;
   float *pfPixels;
   float fAvg;
   float fSum;

	vPrintMsg(4,"GetImageSizesForPrj\n");
	pActImage = imgio_openimage(pchActImageName, 'o', &iXdim, &iYdim, &iZdim);
   if (iXdim != iYdim)
       vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vGetImageSizesForPrj",
         "X and Y dimension of activity image should be equal : x-dimension = %d, Y-dimension = %d", 
         iXdim, iYdim);
	psParms->NumPixels = iXdim;

	iSliceStart  = iGetIntParm("slice_start", &bFound, 0);
	iSliceEnd = iGetIntParm("slice_end", &bFound, iZdim-1); 
	iNumSlices=psParms->NumSlices=iSliceEnd - iSliceStart + 1;
	if (iNumSlices <= 0)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "vGetImageSizesForPrj",
			"Number of slices must be > 0: start=%d, end=%d, num=%d\n",
			iSliceStart, iSliceEnd, psParms->NumSlices);

	vPrintMsg(7,"NumPix=%d, NumSlices=%d \n",
			psParms->NumPixels, psParms->NumSlices);

	imgio_closeimage(pActImage);
}


/**
	 @brief Reads the projection data from the image stored in pImage and puts
	 it into the modified output projection image.
	 
	 Only the requested slices (y values) are put in the image. In the x
	 direction the data are shifted so the COR is in the center of the
	 array. The input image has dimension iNumBins x iNumInSlices x
	 iNumAngles and the output image has dimensions iOutNumBins x
	 iOutNumSlices x iNumAngles.  It is assumed that iStartSlice and
	 iEndSlice have already been checked against the size of the input image.

	 @param *pImage       - image file containing  projection pixels.
	 @param iInNumBins    - number of projection bins (xdim) in input image. 
	 @param iOutNumBins   - number of bins in output modified projections.
	 @param iInNumSlices  - number of slices (ydim) in input image.
	 @param iNumAngles    - number of angles (zdim) in nput and output image.
	 @param iStartSlice   - starting slice (y value) to put at top of output img.
	 @param iOutNumSlices - number of slices (y values) in output image.
	 @param fInBinWidth   - bin width in input prj image file.
	 @param fOutBinWidth  - bin width in the ouput "modified" projections.
	 @param fScaleFac     - factor to scale projections with (ignored if 0 or 1).
	 @param *psViews      - structure containing information about cor for each.
	                        view that is needed to shift the raw projections into
	                        the output matrix.

	 @return ptr to projection data.
*/

static float *pfReadPrjPix(IMAGE *pImage, int iInNumBins, int iOutNumBins,
	int iInNumSlices, int iNumAngles, int iStartSlice, int iOutNumSlices,
	float fInBinWidth, float fOutBinWidth, float fScaleFac, PrjView_t
	*psViews)
{
	float *pfPrjPixels; /* entire set of projection data needed to reconstruct
												 desired slices */
	float *pfReadPix; /* raw projection data for one 2d view */
	int iAngle;
	float fPrjSum=0.0;


	vPrintMsg(6,"ReadPrjPix\n");
	vPrintMsg(7,"            InBins=%d, OutBins=%d, InSlices=%d, OutSlices=%d\n",
		iInNumBins, iOutNumBins, iInNumSlices, iOutNumSlices);
	vPrintMsg(7,"            NumAngles=%d, StartSlice=%d, normfac=%.3g\n",
		iNumAngles, iStartSlice,fScaleFac);
	pfPrjPixels=vector(iOutNumSlices*iOutNumBins*iNumAngles, 
			"ReadPrjImage:PrjPixels");
	set_float(pfPrjPixels, iOutNumSlices*iOutNumBins*iNumAngles, 0.0);
	pfReadPix=vector(iInNumBins*iInNumSlices*iNumAngles,
		"ReadPrjImage:ReadPix");
	for(iAngle=0; iAngle < iNumAngles; ++iAngle){
		/* read all pixels in one projection view */
		imgio_readslices(pImage,iAngle,iAngle,pfReadPix);
		fPrjSum += sum_float(pfReadPix+iStartSlice*iInNumBins, 
							iInNumBins*iOutNumSlices);
		//fprintf(stderr,"%d prjsum=%f\n",iAngle,fPrjSum);
		/* shift only the slices we need into the right place in PrjPixels */
		vMeasToModPrj(iInNumBins, iOutNumBins, iOutNumSlices,
						psViews[iAngle].Left,fInBinWidth,fOutBinWidth,
						pfReadPix+iStartSlice*iInNumBins,
						pfPrjPixels+iAngle*iOutNumSlices*iOutNumBins);
	}
	if (fScaleFac != 0.0 && fScaleFac != 1.0)
		scale_float(pfPrjPixels, iOutNumSlices*iOutNumBins*iNumAngles, 
							fScaleFac);
	vPrintMsg(8,"prjsum=%.2f, modprjsum=%.2f\n",fPrjSum,
		sum_float(pfPrjPixels,iOutNumSlices*iOutNumBins*iNumAngles)); 
	IrlFree(pfReadPix);
	return pfPrjPixels;
}


/**
	 @brief Reads from image and gets the projection slices needed for
	        reconstruction.

	 @note The output image is NOT in "reordered" (y varying slowest) format
	       unless REORDER_PIXELS is defined during compilation. For now the 
			 reordering is performed in irlGenprj()
	 
	 @param *pImage     pointer to input image.
	 @param iStartSlice start slice in input image. 
	 @param iNumSlices  number of slices in output image.
	 @param iNumPixels  number of pixels (xdim and ydim) in input and 
	                    output image.
	 
	 @return a pointer to these pixels.
*/

static float *pfReadIrlImage(IMAGE *pImage, int iStartSlice,
	int iNumSlices, int iNumPixels)
{
	float *pfPixels;
	float *pfReadSlice;
	int iSlice;
	int iX, iY;

	pfPixels = vector(iNumSlices*iNumPixels*iNumPixels, 
		"ReadIrlImage:Pixels");
	
#ifndef REORDER_PIXELS
	imgio_readslices(pImage, iStartSlice, iStartSlice+iNumSlices-1, 
		pfPixels);
#else
	pfReadSlice = vector(iNumPixels*iNumPixels,"ReadIrlImage:ReadSlice");
	for(iSlice=0; iSlice < iNumSlices; ++iSlice){
		imgio_readslices(pImage, iStartSlice+iSlice, iStartSlice+iSlice, 
			pfReadSlice);
		for(iY=0; iY<iNumPixels; ++iY){
			for(iX=0; iX<iNumPixels; ++iX){
				/*store pixels in reordered format*/
				pfPixels[iX + iNumPixels*(iSlice + iY*iNumSlices)]=
					pfReadSlice[iX + iY*iNumPixels];
			}
		}
	}
	IrlFree(pfReadSlice);
#endif
	return(pfPixels);
}

/**
	 @brief Reads the attenunation map.

	 @note The output image is NOT in "reordered" (y varying slowest) format
	       unless REORDER_PIXELS is defined during compilation. For now the
          reordering is performed in irlGenprj()

	 @param *pchAtnMapName - Attenuation data file name.
	 @param *psParms       - Parameters struct.

	 @return ptr to attenuation map.
*/

float *pfGetAtnMap(char *pchAtnMapName, IrlParms_t *psParms)
{
	int iStartSlice, iEndSlice;
	int bFound;
	int iXdim, iYdim, iZdim;
	int i;
	float *pfAtnPixels;
	IMAGE *pAtnImage;

	vPrintMsg(4,"GetAtnMap\n");
	vPrintMsg(6,"atn map=%s\n",pchAtnMapName);
	iStartSlice=iGetIntParm("atn_slice_start",&bFound,0);
	if (iStartSlice < 0)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetAtnMap",
			"Start Slice for Atn Map must be >= 0, not %d",iStartSlice);

	pAtnImage=imgio_openimage(pchAtnMapName,'o',&iXdim,&iYdim,&iZdim);
	if (iYdim != iXdim || iXdim != psParms->NumPixels){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetAtnMap",
			"Atn Map must be same size as activity image");
	}
	iEndSlice = iStartSlice + psParms->NumSlices - 1;
	if (iEndSlice >= iZdim)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetAtnMap",
			"Not enough slices in atn map. First slice=%d, last slice needed: %d",
			iStartSlice, iEndSlice);
	pfAtnPixels = pfReadIrlImage(pAtnImage,iStartSlice, psParms->NumSlices,
											psParms->NumPixels);

#ifdef DEBUG
	writeimage("atnimage.im", psParms->NumPixels, psParms->NumSlices, 
							psParms->NumPixels, pfAtnPixels);
#endif
	imgio_closeimage(pAtnImage);
	return(pfAtnPixels);
}

/**
	 @brief Determines whether there is an additive scatter estimate and if
	 so, reads it in and scale it.
	 
	 @param *psParms
	 @param *psViews

	 @return ptr to scaled scatter estimate.
*/

float *pfGetScatterEstimate(IrlParms_t *psParms, PrjView_t *psViews)
{
	float *pfScatterEstimate=NULL;
	char *pchScatterEstimateFile;
	float fScatFac;
	int iStartSlice, iEndSlice, iNumSlices;
	int iAngle, iNumAngles, iNumBins, iNumPix;
	int bFound;
	int iXdim, iYdim, iZdim;
	IMAGE *pScatImage;

	vPrintMsg(4,"GetScatterEstimate");
	iNumSlices = psParms->NumSlices;
	iNumBins = psParms->NumPixels;
	iNumAngles=psParms->NumViews;

	pchScatterEstimateFile=pchGetStrParm("scat_est_file",&bFound, "");
	if (!bFound){
		vPrintMsg(6,"No scatter estimate file will be used\n");
		return NULL;
	}
	fScatFac = psParms->fScatEstFac=dGetDblParm("scat_est_fac",&bFound,1.0);
	iStartSlice=iGetIntParm("scat_startslice",&bFound, 
						iGetIntParm("slice_start", &bFound, 0));
	iEndSlice = iStartSlice + iNumSlices - 1;
	pScatImage = imgio_openimage(pchScatterEstimateFile, 
							'o', &iXdim, &iYdim, &iZdim);
	if (iXdim != iNumBins || iZdim != iNumAngles)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetScatterEstimate",
			"%s\n%s%s\n",
			"Number of bins or number of angles in scat est file not correct",
			"    Scat est file: " ,pchScatterEstimateFile);
	if (iStartSlice < 0 || iStartSlice > iEndSlice || iEndSlice >= iYdim)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetScatterEstimate",
			"Start & end values for scatter slice are illegal or inconsistent\n");
	/* get the projection pixels, treat them just as we do the projection
		data in terms of shifting*/
	pfScatterEstimate = pfReadPrjPix(pScatImage, iNumBins, iNumBins, iYdim,
		iNumAngles, iStartSlice, iNumSlices, psParms->BinWidth, 
		psParms->BinWidth, fScatFac, psViews);
	imgio_closeimage(pScatImage);
	return (pfScatterEstimate);
}


/**
	 @brief Reads the activity image and extracts the requisite slices.

	 @note The output image is NOT in "reordered" (y varying slowest) format
	       unless REORDER_PIXELS is defined during compilation. For now the
          reordering is performed in irlGenprj()

	 @param pchActImageName
	 @param psParms

	 @return ptr to the activity image.
*/

float *pfGetActImage(char *pchActImageName, IrlParms_t *psParms)
{
	int iStartSlice, iEndSlice;
	int bFound;
	int iXdim, iYdim, iZdim;
	int iX, iY, iZ,i;
	float *pfPixels;
	float fAvg;
	float fSum;
	IMAGE *pActImage;

	vPrintMsg(4,"GetActImage\n");
	iStartSlice=iGetIntParm("slice_start", &bFound, 0);
	if (iStartSlice < 0)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetActImage",
			"Start Slice for Act must be >= 0, not %d",iStartSlice);

	pActImage=imgio_openimage(pchActImageName,'o',&iXdim,&iYdim,&iZdim);
	if (iYdim != iXdim || iXdim != psParms->NumPixels){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetActImage",
			"Activity image must be square (should have been cauth earlier)");
	}
	iEndSlice = iStartSlice + psParms->NumSlices - 1;
	if (iEndSlice >= iZdim)
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetActImage",
		"Not enough slices in activity image. First slice=%d, last slice needed: %d",
			iStartSlice, iEndSlice);
	pfPixels = pfReadIrlImage(pActImage,iStartSlice, psParms->NumSlices,
											psParms->NumPixels);
	imgio_closeimage(pActImage);
	return(pfPixels);
}

/**
	 @brief Reads in the image that is used to compute orbit.

	 @warning The Implementation appears to be INCOMPLETE.

	 @return ptr to the image.
*/

float *pfGetOrbitImage(char *pchFname, IrlParms_t *psParms)
{
	IMAGE *pInImage;
	float *pfPixels;
	float *pfSlice;
	int iXdim, iYdim, iZdim;
	int iImageSize;
	int iSlice, iStartSlice, iEndSlice;
	int i;
	int bFound;

	pInImage=imgio_openimage(pchFname,'o',&iXdim,&iYdim,&iZdim);

	if (iYdim != iXdim || iXdim != psParms->NumPixels){
		vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "ComputeOrbitFromFile",
			"Image used to compute orbit must be square ");
	}
	iImageSize=psParms->NumPixels * psParms->NumPixels, 
	pfPixels = vector(iImageSize, "GetOrbitImage: Pixels");
	if (iZdim == 1){
		//Image is already 2d
		imgio_readslices(pInImage,0,0,pfPixels);
	}else{
		iStartSlice=iGetIntParm("slice_start",&bFound,0);
		iEndSlice = iStartSlice + psParms->NumSlices - 1;
		if (iEndSlice >= iZdim)
			vErrorHandler(ECLASS_FATAL, ETYPE_ILLEGAL_VALUE, "GetOrbitImage",
		"Not enough slices in orbit image. First slice=%d, last slice needed: %d",
			iStartSlice, iEndSlice);
		pfSlice = vector(iImageSize, "GetOrbitImage: Slice");
		set_float(pfPixels,iImageSize,0.0);
		for(iSlice=iStartSlice; i<= iEndSlice; ++i){
			imgio_readslices(pInImage,iSlice,iSlice,pfSlice);
			for(i=0; i<iImageSize; ++i){
				pfPixels[i] += pfSlice[i];
			}
		}
	}
	imgio_closeimage(pInImage);
}
