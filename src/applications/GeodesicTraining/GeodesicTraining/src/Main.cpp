#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaUtilities.h"

#include <iostream>
#include <vector>
#include <unordered_map>

#include "GeodesicTrainingSegmentation.h"

const std::string DEFAULT_MODE = "reversegeotrain";

// Variables need to be global for easier template initialization

cbica::Logging logger;
std::string    loggerFile;
bool           loggerRequested = false;

bool        timerEnabled = false, closeWindow = true, includeDateTime = true, maxThreads = false, loiSet = false, subsample = true, 
            balancedSubsample = GeodesicTrainingSegmentation::DEFAULT_BALANCED_SUBSAMPLE;
short       imageDimensions = 3;
float       threshold = GeodesicTrainingSegmentation::DEFAULT_THRESHOLD, 
            inputImagesToAgdMapsRatio = GeodesicTrainingSegmentation::DEFAULT_INPUT_IMAGES_TO_AGD_MAPS_RATIO;
int         tempPosition, numberOfThreads = 16, maxSamplesForSubsample = GeodesicTrainingSegmentation::DEFAULT_MAX_SAMPLES_SVM_SUBSAMPLE,
            labelOfInterest = GeodesicTrainingSegmentation::DEFAULT_LABEL_OF_INTEREST,
            labelTC = GeodesicTrainingSegmentation::DEFAULT_LABEL_TC, labelET = GeodesicTrainingSegmentation::DEFAULT_LABEL_ET,
            labelED = GeodesicTrainingSegmentation::DEFAULT_LABEL_ED, labelHT = GeodesicTrainingSegmentation::DEFAULT_LABEL_HT;
std::string labelsPath = "", mode = DEFAULT_MODE, outputDir = "./", configFilePath = "", rfConfigFilePath = "",
            datasetName = "", tag = "", pixelType = "float", groundTruthPath = "",
            flairPath = "", t1Path = "", t1cePath = "", t2Path = "";

std::vector< std::string >     input, inputModels;
std::vector< int >             groundTruthSkip;
std::vector< double >          importanceValues;
std::unordered_map< int, int > changeLabelsMap;

// Functions

void Run();

template <typename PixelType, unsigned int Dimensions> 
void RunGeodesicTraining();

int main(int argc, char *argv[])
{
	cbica::CmdParser parser = cbica::CmdParser(argc, argv, "GeodesicTraining");

	// Optional Parameters

	std::string modeMsg = "Available modes are: \n\n";
	modeMsg += "                             reversegeotrain -> [DEFAULT MODE] use AGD and SVMs to produce labels\n";
	modeMsg += "                             svmlabels       -> use SVMs to produce labels\n";
	modeMsg += "           [2-class]         svmpseudo       -> use SVMs to produce pseudoprobability maps\n";
	modeMsg += "           [1-class]         agd             -> run AGD to produce agd maps\n";
	modeMsg += "           [2-class]         geotrain        -> use SVMs and AGD to produce agd maps and threshold them\n";
	modeMsg += "           [2-class]         geotrainfull    -> AGD ---> SVMs ---> AGD ---> threshold\n";
	modeMsg += "                             rf              -> Random Forest\n";
	modeMsg += "                             rfauto          -> Random Forest with automatic parameter tuning\n";
	modeMsg += "                             agdrf           -> Random Forest with AGD input images\n";
	modeMsg += "                             agdrfauto       -> Random Forest with AGD input images and automatic parameter tuning\n";
	modeMsg += "           [util]            segment         -> segment an area from an image according to a labels image\n";
	modeMsg += "           [util]            labelthres      -> convert AGD maps to labels images\n";
	modeMsg += "           [util]            checkaccuracy   -> Check the accuracy of a segmentation in relation to a ground truth one\n";
	modeMsg += "           [util]            changelabels    -> change the labels of a labels image\n";
	modeMsg += "           [util]            generateconfig  -> generate a configuration file from pretrained svm models\n";
	parser.addOptionalParameter("m", "mode", cbica::Parameter::STRING, "text", modeMsg);

	parser.addOptionalParameter("i", "input", cbica::Parameter::STRING, ".nii.gz", "List of full paths to the input files, separated by comma.");
	parser.addOptionalParameter("L", "Logger", cbica::Parameter::STRING, "log file with user write access", "Full path to log file to store console output.",
		"By default, only console output is generated.");
	parser.addOptionalParameter("l", "labels", cbica::Parameter::STRING, ".nii.gz", "Full path to the input labels file.");
	parser.addOptionalParameter("o", "output", cbica::Parameter::STRING, ".nii.gz", "Directory for output.");
	parser.addOptionalParameter("c", "configuration", cbica::Parameter::STRING, ".yaml", "Full path to SVM Configuration file.");
	parser.addOptionalParameter("rfc", "rfconfiguration", cbica::Parameter::STRING, ".config", "Full path to RF Configuration file.");
	parser.addOptionalParameter("d", "datasetname", cbica::Parameter::STRING, "text", "The dataset name", "Used for better output organization.",
		"Folders will be created inside the output dir.");
	parser.addOptionalParameter("tg", "tag", cbica::Parameter::STRING, "text", "A tag for this execution", "Used for better output organization.",
		"Tag will be used only if dataset name is set."); 
	parser.addOptionalParameter("nd", "nodatetime", cbica::Parameter::STRING, "-", "Don't include datetime in the directory name", 
		"Used for better output organization.", "Datetime will still be used if datasetname is provided but no tag.");
	parser.addOptionalParameter("n", "noclose", cbica::Parameter::STRING, "-", "Option not to close the window after the program has finished.");
	parser.addOptionalParameter("r", "reportseconds", cbica::Parameter::STRING, "-", "Enable time report.",
		"Creates time_report.txt file in output directory.");
	parser.addOptionalParameter("t", "threshold", cbica::Parameter::STRING, "float", "Threshold value");
	parser.addOptionalParameter("cl", "changelabelsmap", cbica::Parameter::STRING, "list of tuples", "Importance for each model separated by comma.",
		"Example: Change label 1 -> 2, and labels 3 and 4 -> 5",
		"the changelabelsmap value should be \"1/2,3/5,4/5\"");
	parser.addOptionalParameter("im", "inputmodels", cbica::Parameter::STRING, ".xml", "List of input models separated by comma.",
		"For mode: \"generateconfig\"");
	parser.addOptionalParameter("iv", "importancevalue", cbica::Parameter::STRING, "double vector", "Importance for each model separated by comma.",
		"For mode: \"generateconfig\"");
	parser.addOptionalParameter("g", "groundtruth", cbica::Parameter::STRING, ".nii.gz", "Ground truth image (for checking accuracy).");
	parser.addOptionalParameter("gs", "groundtruthskip", cbica::Parameter::STRING, "int", "List of labels to skip when checking accuracy ",
		"in relation to ground truth image. Default is 0.");
	parser.addOptionalParameter("j", "threads", cbica::Parameter::STRING, "int", "Number of threads that will be used.", 
		"Give \"max\" for infinite.", "Default is 16.");
	parser.addOptionalParameter("flair", "flair", cbica::Parameter::STRING, ".nii.gz", "Flair image (For MRI images only)",
		"Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images.");
	parser.addOptionalParameter("t1", "t1", cbica::Parameter::STRING, ".nii.gz", "T1 image (For MRI images only)",
		"Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images.");
	parser.addOptionalParameter("t1ce", "t1ce", cbica::Parameter::STRING, ".nii.gz", "T1ce image (For MRI images only)",
		"Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images.");
	parser.addOptionalParameter("t2", "t2", cbica::Parameter::STRING, ".nii.gz", "T2 image (For MRI images only)",
		"Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images.");
	parser.addOptionalParameter("ltc", "labeltc", cbica::Parameter::STRING, "int", "Label for tumor core (For MRI images only)");
	parser.addOptionalParameter("let", "labelet", cbica::Parameter::STRING, "int", "Label for enhanced tumor (For MRI images only)");
	parser.addOptionalParameter("led", "labeled", cbica::Parameter::STRING, "int", "Label for edema (For MRI images only)");
	parser.addOptionalParameter("lht", "labelht", cbica::Parameter::STRING, "int", "Label for healthy tissue (For MRI images only)");
	parser.addOptionalParameter("loi", "labelofinterest", cbica::Parameter::STRING, "int", "Label of interest.",
		"Default is 1.", "Any other non-zero label can be used for areas that are not of interest.");
	parser.addOptionalParameter("ns", "nosubsample", cbica::Parameter::STRING, "-", "Option to not subsample if the number of samples is high (SVM)");
	//parser.addOptionalParameter("ir", "imagestoagdratio", cbica::Parameter::STRING, "int", "Input images to AGD maps ratio (default is 6)",
	//	"The bigger the value, the bigger the effect of", "the input images compared to the AGD maps");
	parser.addOptionalParameter("ms", "maxsamples", cbica::Parameter::STRING, "int", "Max samples to use (SVM Subsampling)");
	parser.addOptionalParameter("nb", "nobalancesamples", cbica::Parameter::STRING, "int", "Don't balance the subsampling (SVM Subsampling)");
	//parser.addOptionalParameter("pt", "pixeltype", cbica::Parameter::STRING, "string", "Input image(s) pixel type (int or float) [default is float]");
	parser.addOptionalParameter("id", "imagedimensions", cbica::Parameter::STRING, "int", "Input image(s) dimensions [only 3D supported for now]");

	// Parameters parsing

	if ((argc < 1) || (parser.compareParameter("u", tempPosition))) {
		parser.echoUsage();
		return EXIT_SUCCESS;
	}
	if (parser.compareParameter("h", tempPosition)) {
		parser.echoHelp();
		return EXIT_SUCCESS;
	}
	if (parser.compareParameter("v", tempPosition)) {
		parser.echoVersion();
		return EXIT_SUCCESS;
	}
	if (parser.compareParameter("L", tempPosition)) {
		loggerFile = argv[tempPosition + 1];
		loggerRequested = true;
		logger.UseNewFile(loggerFile);
	}
	if (parser.compareParameter("m", tempPosition)) {
		// Mode
		mode = argv[tempPosition + 1];
	}
	if (parser.compareParameter("i", tempPosition)) {
		// Input files
		std::string inputPathsRaw = argv[tempPosition + 1];
		input = cbica::stringSplit(inputPathsRaw, ",");
	}
	if (parser.compareParameter("l", tempPosition)) {
		// Labels file
		labelsPath = argv[tempPosition + 1];
	}
	if (parser.compareParameter("o", tempPosition)) {
		// Output directory
		outputDir = argv[tempPosition + 1];
	}
	if (parser.compareParameter("c", tempPosition)) {
		// Config file
		configFilePath = argv[tempPosition + 1];
	}
	if (parser.compareParameter("rfc", tempPosition)) {
		// Config file
		rfConfigFilePath = argv[tempPosition + 1];
	}
	if (parser.compareParameter("d", tempPosition)) {
		// Dataset name
		datasetName = argv[tempPosition + 1];
	}
	if (parser.compareParameter("tg", tempPosition)) {
		// Tag
		tag = argv[tempPosition + 1];
	}
	if (parser.compareParameter("nd", tempPosition)) {
		// No datetime
		includeDateTime = false;
	}
	if (parser.compareParameter("n", tempPosition)) {
		// Keep window open after execution
		closeWindow = false;
	}
	if (parser.compareParameter("r", tempPosition)) {
		// Create time report
		timerEnabled = true;
	}
	if (parser.compareParameter("t", tempPosition)) {
		// Threshold
		threshold = std::stof(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("im", tempPosition)) {
		// List of input models (For mode "generateconfig")
		std::string imsRaw = argv[tempPosition + 1];
		inputModels = cbica::stringSplit(imsRaw, ",");
	}
	if (parser.compareParameter("iv", tempPosition)) {
		// List of importance values (For mode "generateconfig")
		std::string ivsRaw = argv[tempPosition + 1];
		std::vector< std::string > ivsStrings = cbica::stringSplit(ivsRaw, ",");

		for (auto ivString : ivsStrings) {
			importanceValues.push_back(std::stod(ivString));
		}
	}
	if (parser.compareParameter("g", tempPosition)) {
		// Ground truth: image path
		groundTruthPath = argv[tempPosition + 1];
	}
	if (parser.compareParameter("gs", tempPosition)) {
		// Ground truth: label to skip when comparing
		std::string gsRaw = argv[tempPosition + 1];
		auto gsStrings = cbica::stringSplit(gsRaw, ",");

		for (auto gsString : gsStrings) {
			groundTruthSkip.push_back(stoi(gsString));
		}
	}
	if (parser.compareParameter("j", tempPosition)) {
		// Number of threads
		std::string n = argv[tempPosition + 1];
		if (n == "max") {
			maxThreads = true;
		}
		else if (n != "") {
			numberOfThreads = std::stoi(n);
		}
	}
	if (parser.compareParameter("flair", tempPosition)) {
		// MRI
		flairPath = argv[tempPosition + 1];
	}
	if (parser.compareParameter("t1", tempPosition)) {
		// MRI
		t1Path = argv[tempPosition + 1];
	}
	if (parser.compareParameter("t1ce", tempPosition)) {
		// MRI
		t1cePath = argv[tempPosition + 1];
	}
	if (parser.compareParameter("t2", tempPosition)) {
		// MRI
		t2Path = argv[tempPosition + 1];
	}
	if (parser.compareParameter("ltc", tempPosition)) {
		// MRI
		labelTC = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("let", tempPosition)) {
		// MRI
		labelET = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("led", tempPosition)) {
		// MRI
		labelED = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("lht", tempPosition)) {
		// MRI
		labelHT = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("loi", tempPosition)) {
		// For 2-class (and alternative to ltc/let/led/lht)
		loiSet = true;
		labelOfInterest = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("ns", tempPosition)) {
		// No subsample
		subsample = false;
	}
	//if (parser.compareParameter("ir", tempPosition)) {
	//	// Input images to AGD maps ratio
	//	inputImagesToAgdMapsRatio = std::stof(argv[tempPosition + 1]);
	//}
	if (parser.compareParameter("ms", tempPosition)) {
		// Max samples (if there are more samples than max, only max will be kept)
		maxSamplesForSubsample = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("nb", tempPosition)) {
		// Don't balance the subsampling. Train test split will be completely random
		balancedSubsample = false;
	}
	//if (parser.compareParameter("pt", tempPosition)) {
	//	// Input image(s) pixel type (short, float or int)
	//	pixelType = argv[tempPosition + 1];
	//}
	if (parser.compareParameter("id", tempPosition)) {
		// Input image(s) dimensions
		imageDimensions = std::stoi(argv[tempPosition + 1]);
	}
	if (parser.compareParameter("cl", tempPosition)) {
		// Change labels list (The label before the '/' will be changed with the label after the '/')
		std::string clsRaw = argv[tempPosition + 1];

		if (clsRaw != "") {
			std::vector< std::string > clTuplesRaw = cbica::stringSplit(clsRaw, ",");
			for (auto clTupleRaw : clTuplesRaw) {
				auto clTuple = cbica::stringSplit(clTupleRaw, "/");
				changeLabelsMap[std::stoi(clTuple[0])] = std::stoi(clTuple[1]);
			}
		}
	}
	else {
		// Default change labels if no -cl parameter is provided and at least one of FLAIR/T1/T1CE/T2 is set.
		if (flairPath != "" || t1Path != "" || t1cePath != "" || t2Path != "") {
			changeLabelsMap[labelHT] = 0;
		}
	}

	// Check pixel type and dimensions

	if (pixelType != "float") {
		// This is kept because other types might be supported in the future
		std::cerr << "Only float pixel type images are supported\n";
		exit(EXIT_FAILURE);
	}

	if (imageDimensions != 2 && imageDimensions != 3) {
		std::cerr << "Only 2D or 3D Images supported\n";
		exit(EXIT_FAILURE);
	}

	// Calling the geodesic training algorithm
	Run();

	// Finishing
	if (!closeWindow) {
		std::cout << "\nProgram finished. Press enter to continue...";
		std::cin.get();
	}
	return 0;
}

/** This function is needed, otherwise code will be duplicated for each combination of pixel type and image dimensions */
void Run()
{
	if (pixelType == "float" && imageDimensions == 2) {
		RunGeodesicTraining<float, 2>();
	}
	else if (pixelType == "float" && imageDimensions == 3) {
		RunGeodesicTraining<float, 3>();
	}
}

/** Passes the parameters and executes the algorithm */
template <typename PixelType, unsigned int Dimensions>
void RunGeodesicTraining()
{
	GeodesicTrainingSegmentation::Coordinator<PixelType, Dimensions> geodesicTraining;

	// For mode

	if (mode == "svmpseudo") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::SVM_PSEUDO);
	}
	else if (mode == "svmlabels") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::SVM_LABELS);
	}
	else if (mode == "agd") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::AGD);
	}
	else if (mode == "labelsthres") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::LABELS_THRESHOLD);
	}
	else if (mode == "segment") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::SEGMENT);
	}
	else if (mode == "reversegeotrain") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::REVERSE_GEOTRAIN);
	}
	else if (mode == "geotrain") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::GEOTRAIN);
	}
	else if (mode == "geotrainfull") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::GEOTRAIN_FULL);
	}
	else if (mode == "rf") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::RF);
	}
	else if (mode == "agdrf") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::AGD_RF);
	}
	else if (mode == "rfauto") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::RF_AUTO);
	}
	else if (mode == "agdrfauto") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::AGD_RF_AUTO);
	}
	else if (mode == "changelabels") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::CHANGE_LABELS);
	}
	else if (mode == "generateconfig") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::GENERATE_CONFIG);
	}
	else if (mode == "checkaccuracy") {
		geodesicTraining.SetMode(GeodesicTrainingSegmentation::MODE::CHECK_ACCURACY);
	}
	else {
		std::cerr << "Error: Invalid mode.\n";
		exit(EXIT_FAILURE);
	}

	// For special MRI modalities

	geodesicTraining.SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI::FLAIR, flairPath);
	geodesicTraining.SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI::T1, t1Path);
	geodesicTraining.SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI::T1CE, t1cePath);
	geodesicTraining.SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI::T2, t2Path);
	geodesicTraining.SetTumorCoreLabelMRI(labelTC);
	geodesicTraining.SetEnhancedTumorLabelMRI(labelET);
	geodesicTraining.SetEdemaLabelMRI(labelED);
	geodesicTraining.SetHealthyTissueLabelMRI(labelHT);

	// For other parameters

	geodesicTraining.SetVerbose(true);
	geodesicTraining.SetInputImages(input);
	geodesicTraining.SetLabels(labelsPath);
	geodesicTraining.SetConfigFile(configFilePath);     // For modes that use SVMs
	geodesicTraining.SetRfConfigFile(rfConfigFilePath); // For modes that use RFs
	geodesicTraining.SetGroundTruth(groundTruthPath, groundTruthSkip);
	geodesicTraining.SetChangeLabelsMap(changeLabelsMap);
	geodesicTraining.SetThreshold(threshold);
	//geodesicTraining.SetInputImageToAgdMapsRatio(inputImagesToAgdMapsRatio);
	geodesicTraining.SetOutputPath(outputDir, datasetName, tag, includeDateTime);
	geodesicTraining.SetSaveAll(true);
	geodesicTraining.SetTimerEnabled(timerEnabled);
	geodesicTraining.SetNumberOfThreads(numberOfThreads);
	geodesicTraining.SetNumberOfThreadsMax(maxThreads);
	geodesicTraining.SetSubsampling(subsample, maxSamplesForSubsample);
	geodesicTraining.SetPretrainedModelsPaths(inputModels); // For mode "generateconfig"
	geodesicTraining.SetImportanceValues(importanceValues); // For mode "generateconfig"

	if (loiSet) {
		geodesicTraining.SetLabelOfInterest(labelOfInterest);
	}

	// Executing
	geodesicTraining.Execute();
}