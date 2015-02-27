#define M_PI  3.141592653

/* genprj.c */
void vSetCallBackData(int iNumPixels, int iNumSlices, int iNumAngles, char *pchOutBase, PrjView_t *pView, IrlParms_t *psParm);
void vPrjViewCallback(int iAngle, float *pfCalcPrj);
char *pUsageMsg(void);
int iSetupFromCmdLine(int iArgc, char **ppchArgv, IrlParms_t *psIrlParms, Options_t *psOptions, PrjView_t **ppsPrjViews, char **ppchSrfKrnlFile, char **ppchDrfTabFile, char **ppchLogFile, char **ppchMsgFile, char **ppchOutBase, float **ppfPrjImage, float **ppfAtnMap, float **ppfReconImage, float **ppfScatterEstimate);
int main(int iArgc, char **ppchArgv);

/* GenprjSetup.c */
char *pchGetNormBase(char *pchBase);
void vGetPrjParms(IrlParms_t *psParms, Options_t *psOptions);
int iParseModelString(char *pchModelStr, char *pchName);
void vGetEffectsToModel(int *pbModelAtn, int *pbModelDrf, int *pbModelSrf);
PrjView_t *psSetupPrjViews(IrlParms_t *psParms);
char *pchGetSrfKrnlFname(void);
char *pchGetDrfTabFname(void);

/* GetImages.c */
void vGetImageSizesForPrj(char *pchPrjImageName, IrlParms_t *psParms);
float *pfGetAtnMap(char *pchAtnMapName, IrlParms_t *psParms);
float *pfGetScatterEstimate(IrlParms_t *psParms, PrjView_t *psViews);
float *pfGetActImage(char *pchActImageName, IrlParms_t *psParms);
float *pfGetOrbitImage(char *pchFname, IrlParms_t *psParms);

/* MeasToModPrj.c */
void vMeasToModPrj(int nBins, int nRotPixs, int nSlices, float Left, float BinWidth, float PixelWidth, float *RawPrjData, float *ModPrjData);
void Interp_bck( float *InPrjData, float *OutPrjData, IrlParms_t *psParms, PrjView_t psView, float *Sum);
