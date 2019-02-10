#ifndef H_CBICA_RF_PREPARE_TRAIN_DATA
#define H_CBICA_RF_PREPARE_TRAIN_DATA

#include <opencv2/ml.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <vector>

namespace RFSuiteTrainData
{
	/**
	Creates TrainData where random x rows will be the training set (and random (data.rows-x) will be the testing set)
	[NOT USED]
	@param data the input matrix
	@param responses the input input responses
	@param ntrain_samples the x
	@return pointer to the TrainData object
	*/
	cv::Ptr<cv::ml::TrainData>
	PrepareTrainData(const cv::Mat &data, const cv::Mat &responses, int ntrain_samples);

	/**
	Shuffles rows of two "parallel" cv::Mat matrices
	[NOT USED]
	@param matrix the first matrix
	@param resRandMatrix the result randomized first matrix
	@param responses the second matrix
	@param resRandResponses the result randomized second matrix
	*/
	void 
	shuffleDataAndResponses(const cv::Mat &matrix, cv::Mat &resRandMatrix, const cv::Mat &responses, cv::Mat &resRandResponses);

	/**
	Shuffles rows of a cv::Mat
	[NOT USED]
	@param matrix the input matrix
	@return the shuffled matrix
	*/
	cv::Mat 
	shuffleRows(const cv::Mat &matrix);
}

#endif // !H_CBICA_RF_PREPARE_TRAIN_DATA