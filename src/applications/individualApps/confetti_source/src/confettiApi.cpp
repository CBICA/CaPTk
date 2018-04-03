
#include "confetti.h"
#include "multinomialModel.h"
#ifdef _WIN32
# define DllExport  extern "C" __declspec( dllexport ) 
#else
#define DllExport extern "C"
#endif

Confetti g_core;
static TXT_CALLBACK_TYPE g_txtCallback;
DllExport void registerCallback( TXT_CALLBACK_TYPE f )
{
	g_txtCallback = f;
	g_core.setListner(f);
}
DllExport bool genConnectivitySig(const char* tdiCsvFile, const char* signatureFileName, const char* fiberFile)
{
	g_core.genConnectivitySig(tdiCsvFile, signatureFileName, fiberFile);
	return true;
}
DllExport bool genAdaptiveCluster(const char* signatureFile, const char* tempRoot, const char* outFile)
{
	MNM obj;
	obj.setListner(g_txtCallback);
	obj.genAdaptiveCluster(signatureFile, tempRoot, outFile);
	return true;
}
DllExport bool extractTracts(const char* outDIR, const char* tempRoot, const char* inputIDFile, const char* inputNamesFile, const char* inputLabelFile, const char* inputBfloatFile)
{
	g_core.extractTract(outDIR, tempRoot, inputIDFile, inputNamesFile, inputLabelFile, inputBfloatFile);
	return true;
}

