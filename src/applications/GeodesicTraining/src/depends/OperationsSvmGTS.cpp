#include "OperationsSvmGTS.h"

bool GeodesicTrainingSegmentation::CreateBalancedSubsample(std::shared_ptr<ParserGTS::Result>& data, 
		std::string& errorMessageIfApplicable,
		std::unordered_map<int, int> labelsCountMap, int maxSamples)
	{
		typedef int LabelsPixelType;

		// Find how many samples to keep for each label
		std::unordered_map<LabelsPixelType, int> maxForEachLabel;
		double totalLabels = 0;

		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			// Initialize
			maxForEachLabel[*it] = 0;
		}

		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			totalLabels += labelsCountMap[*it];
		}

		std::cout << "\n";
		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			maxForEachLabel[*it] = std::lround(maxSamples * labelsCountMap[*it] / totalLabels);
			std::cout << "MAX FOR " << *it << ": " << maxForEachLabel[*it] << "\n";
			if (maxForEachLabel[*it] < 3) {
				errorMessageIfApplicable = std::string("Cannot create a balanced subsample, ") +
					std::string("because one (or more) labels is outweighted by the others a lot. ") +
					std::string("Draw more for the labels that you haven't drawn a lot for, or remove ") +
					std::string("some samples from the ones that are large.");
				return false;
			}
		}

		// Shuffle training matrix and labels (with the same seed)
		std::vector <int> seeds;
		for (int cont = 0; cont < data->trainingMat.rows; cont++) {
			seeds.push_back(cont);
		}

		cv::randShuffle(seeds);

		cv::Mat randTrainingMat, randLabelsMat;
		for (int cont = 0; cont < data->trainingMat.rows; cont++) {
			randTrainingMat.push_back(data->trainingMat.row(seeds[cont]));
			randLabelsMat.push_back(data->labelsMat.row(seeds[cont]));
		}
		
		data->trainingMat = randTrainingMat;
		data->labelsMat = randLabelsMat;

		// // Create a mask for which samples to keep for training
		// cv::Mat sampleIdx = cv::Mat::zeros(data->trainingMat.rows, 1, CV_8U);

		std::unordered_map<LabelsPixelType, int> keptForEachLabel;

		for (std::set<LabelsPixelType>::iterator it = data->differentLabels.begin(); it != data->differentLabels.end(); it++)
		{
			// Initialize
			keptForEachLabel[*it] = 0;
		}

		cv::Mat trainingMatNew, labelsMatNew;

		for (int row_i = 0; row_i < data->trainingMat.rows; row_i++)
		{
			LabelsPixelType label = data->labelsMat.ptr<LabelsPixelType>(row_i)[0];
			
			if (keptForEachLabel[label] <= maxForEachLabel[label]) {
				//sampleIdx.ptr<uchar>(row_i)[0] = 1;
				keptForEachLabel[label] += 1;
				trainingMatNew.push_back(data->trainingMat.row(row_i).clone());
				labelsMatNew.push_back(data->labelsMat.row(row_i).clone());
			}
		}
		data->trainingMat = trainingMatNew;
		data->labelsMat   = labelsMatNew;

		return true;
		//return sampleIdx;
	}