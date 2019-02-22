#include "ConfigParserRF.h"

void ConfigParserRF::Parse(std::string filePath, double & trainingSamplePercentage, int & maxDepth, double & minSampleCountPercentage, int & maxCategories, int & activeVarCount, int & numberOfTrees, cv::Mat &priors)
{
	std::ifstream infile(filePath);

	infile >> trainingSamplePercentage >> maxDepth >> minSampleCountPercentage >> maxCategories >> activeVarCount >> numberOfTrees;

	float val;
	while (infile >> val) {
		if (val == 0) {
			break; // No priors will be used
		}
		priors.push_back(val);
	}
}

void ConfigParserRF::PrintParseResult(double trainingSamplePercentage, int maxDepth, double minSampleCountPercentage, int maxCategories, int activeVarCount, int numberOfTrees, cv::Mat & priors)
{
	std::cout << "\nRF CONFIG:";
	std::cout << "\n\tTRAINING SAMPLE %:  " << trainingSamplePercentage;
	std::cout << "\n\tMAX DEPTH:          " << maxDepth;
	std::cout << "\n\tMIN SAMPLE COUNT %: " << minSampleCountPercentage;
	std::cout << "\n\tMAX CATEGORIES:     " << maxCategories;
	std::cout << "\n\tACTIVE VAR COUNT:   " << activeVarCount;
	std::cout << "\n\tNUMBER OF TREES:    " << numberOfTrees;
	std::cout << "\n\tPRIORS:";

	for (int i = 0; i < priors.rows; i++)
	{
		std::cout << " " << priors.ptr<float>(i)[0];
	}
	std::cout << "\n\n";
}

void ConfigParserRF::PrintParseResultToFile(std::string filePath, double trainingSamplePercentage, int maxDepth, double minSampleCountPercentage, int maxCategories, int activeVarCount, int numberOfTrees, cv::Mat & priors)
{
	std::ofstream rfReportFile;
	rfReportFile.open(filePath, std::ios_base::app); //append file

	rfReportFile << "\nRF CONFIG:";
	rfReportFile << "\n\tTRAINING SAMPLE %:  " << trainingSamplePercentage;
	rfReportFile << "\n\tMAX DEPTH:          " << maxDepth;
	rfReportFile << "\n\tMIN SAMPLE COUNT %: " << minSampleCountPercentage;
	rfReportFile << "\n\tMAX CATEGORIES:     " << maxCategories;
	rfReportFile << "\n\tACTIVE VAR COUNT:   " << activeVarCount;
	rfReportFile << "\n\tNUMBER OF TREES:    " << numberOfTrees;
	rfReportFile << "\n\tPRIORS:";

	for (int i = 0; i < priors.rows; i++)
	{
		rfReportFile << " " << priors.ptr<float>(i)[0];
	}
	rfReportFile << "\n\n";
}