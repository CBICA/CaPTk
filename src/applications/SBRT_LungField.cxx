#include <time.h>
#include <string>
#include "SBRT_LungField.h"

#include "cbicaCmdParser.h"

using namespace std;

int main( int argc,char**argv )
{
	clock_t processStart = clock();

	const int imageDimension = 3;	// image dimension, now support 3D only
	const int inputImageNum = 2;	// number of input images (modality)

	if (imageDimension != 3)
	{
		cerr << " Support 3D image only !" << endl;

		return EXIT_FAILURE;
	}

	cbica::CmdParser parser(argc, argv);
	parser.addRequiredParameter("p", "petImage", cbica::Parameter::STRING, "none", "Absolute path of PET image", "For example /data/.../pet.nii.gz");	
	parser.addRequiredParameter("c", "ctImage", cbica::Parameter::STRING, "none", "Absolute path of CT image", "For example /data/.../ct.nii.gz");
	parser.addRequiredParameter("o", "outputImage", cbica::Parameter::STRING, "none", "Absolute path and basename of output file (without extension)", "For example /output/.../label");
	parser.addOptionalParameter("m", "maskImage", cbica::Parameter::STRING, "none", "Absolute path of mask image", "For example /data/.../mask.nii.gz");
	parser.addOptionalParameter("L", "logFile", cbica::Parameter::STRING, "none", "Absolute path of log file", "For example log_file.txt");
	parser.exampleUsage("SBRT_LungField -p AAAC0_flair_pp_shrunk_75.nii.gz -c AAAC0_flair_pp_shrunk_75.nii.gz -o <output dir>");

	string petName;
	string ctName;
	string oName;
	string maskName;
	string logName;

	if (parser.isPresent("p"))
	{
		parser.getParameterValue("p", petName);
	}
	if (parser.isPresent("c"))
	{
		parser.getParameterValue("c", ctName);
	}
	if (parser.isPresent("o"))
	{
		parser.getParameterValue("o", oName);
	}
	if (parser.isPresent("m"))
	{
		parser.getParameterValue("m", maskName);
	}
	if (parser.isPresent("L"))
	{
		parser.getParameterValue("L", logName);
	}

	vector<string> inputFileName;
	inputFileName.push_back(petName);
	inputFileName.push_back(ctName);
	//for ( int i=0;i<inputImageNum;i++ ) // PET first, CT second
	//{
	//	inputFileName.push_back( argv[i+1] );
	//}

    //string oName = argv[inputImageNum+1];
    
	//float svxSize = atof( argv[inputImageNum+2] );
	//float compactness = atof( argv[inputImageNum+3] );
	//float minSize = atof( argv[inputImageNum+4] );
	//int iter = atoi( argv[inputImageNum+5] );
	//int ptSz = atoi( argv[inputImageNum+6] );
	float svxSize = 512;
	float compactness = 1;
	float minSize = 125;
	int iter = 150;
	int ptSz = 1;

	int maskAvail = 0;
	if (!maskName.empty())
		maskAvail = 1;
	
	//
	SBRT_LungField< float,imageDimension,inputImageNum > lfObject;

	if (!logName.empty())
	{
		lfObject.SetLogger(logName);
	}

	lfObject.SetParameters(maskAvail,svxSize,compactness,minSize,iter,ptSz );
	if (maskAvail==1)
	{
		lfObject.ReadMask(maskName);
	}
	lfObject.Initialization( inputFileName );
	lfObject.DoSupervoxelSegmentation();
	lfObject.GetSvxlFea( inputFileName[1] );
	lfObject.DoKmeans( 5, inputFileName[1] );

  lfObject.WriteLabel(oName);// , maskAvail );

	clock_t processEnd = clock();
	double processDuration = (double)(processEnd - processStart)/CLOCKS_PER_SEC;
	cout << "the whole process costs: " << processDuration << "s" << endl;

	return 0;
}
