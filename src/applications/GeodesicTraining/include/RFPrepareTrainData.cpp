#include "RFPrepareTrainData.h"

cv::Ptr<cv::ml::TrainData> RFSuiteTrainData::PrepareTrainData(const cv::Mat &data, const cv::Mat &responses, int ntrain_samples)
{
	// Mask for which of the first ntrain_samples of the data to be chosen as training data
	cv::Mat sample_idx = cv::Mat::zeros(1, data.rows, CV_8U);
	cv::Mat train_samples = sample_idx.colRange(0, ntrain_samples);
	train_samples.setTo(cv::Scalar::all(1));

	// The types (see definition of TrainData)
	int nvars = data.cols;
	cv::Mat var_type(nvars + 1, 1, CV_8U);
	var_type.setTo(cv::Scalar::all(cv::ml::VAR_ORDERED));
	var_type.at<uchar>(nvars) = cv::ml::VAR_CATEGORICAL;

	// Randomize data with the same seed
	cv::Mat randData, randResp;
	shuffleDataAndResponses(data, randData, responses, randResp);

	return cv::ml::TrainData::create(randData, cv::ml::ROW_SAMPLE, randResp,
		cv::noArray(), sample_idx, cv::noArray(), var_type);
}

void RFSuiteTrainData::shuffleDataAndResponses(const cv::Mat &matrix, cv::Mat &resRandMatrix, const cv::Mat &responses, cv::Mat &resRandResponses)
{
	std::vector <int> seeds;
	for (int cont = 0; cont < matrix.rows; cont++) {
		seeds.push_back(cont);
	}

	cv::randShuffle(seeds);

	for (int cont = 0; cont < matrix.rows; cont++) {
		resRandMatrix.push_back(matrix.row(seeds[cont]));
		resRandResponses.push_back(responses.row(seeds[cont]));
	}
}

cv::Mat RFSuiteTrainData::shuffleRows(const cv::Mat &matrix)
{
	std::vector <int> seeds;
	for (int cont = 0; cont < matrix.rows; cont++)
		seeds.push_back(cont);

	cv::randShuffle(seeds);

	cv::Mat output;
	for (int cont = 0; cont < matrix.rows; cont++)
		output.push_back(matrix.row(seeds[cont]));

	return output;
}