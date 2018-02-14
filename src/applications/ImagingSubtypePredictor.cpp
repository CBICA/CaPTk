#include "ImagingSubtypePredictor.h"


typedef itk::Image< float, 3 > ImageType;
//ImagingSubtypePredictor::~ImagingSubtypePredictor()
//{
//}

VectorDouble ImagingSubtypePredictor::GetStatisticalFeatures(const VectorDouble &intensities)
{
  VectorDouble StatisticalFeatures;

  double temp = 0.0;
  double mean = 0.0;
  double std = 0.0;

  for (unsigned int i = 0; i < intensities.size(); i++)
    temp = temp + intensities[i];
  mean = temp / intensities.size();

  temp = 0;
  for (unsigned int i = 0; i < intensities.size(); i++)
    temp = temp + (intensities[i] - mean)*(intensities[i] - mean);
  std = std::sqrt(temp / (intensities.size() - 1));

  StatisticalFeatures.push_back(mean);
  StatisticalFeatures.push_back(std);

  return StatisticalFeatures;
}

VectorDouble ImagingSubtypePredictor::GetHistogramFeatures(const VectorDouble &intensities)
{
	VectorDouble BinCount;
	for (int i = 0; i <= 4; i++)
		BinCount.push_back(0);

  for (unsigned int i = 0; i < intensities.size(); i++)
  {
	  if (intensities[i] >= 0 && intensities[i] <= 49)
		  BinCount[0] = BinCount[0] + 1;
	  else if (intensities[i] >= 50 && intensities[i] <= 99)
		  BinCount[1] = BinCount[1] + 1;
	  if (intensities[i] >= 100 && intensities[i] <= 149)
		  BinCount[2] = BinCount[2] + 1;
	  if (intensities[i] >= 150 && intensities[i] <= 199)
		  BinCount[3] = BinCount[3] + 1;
	  if (intensities[i] >= 200)
		  BinCount[4] = BinCount[4] + 1;
  }
  return BinCount;
}

VectorDouble ImagingSubtypePredictor::GetVolumetricFeatures(const double &edemaSize,const double &tuSize, const double &neSize, const double &totalSize)
{
  VectorDouble VolumetricFeatures;
  VolumetricFeatures.push_back(edemaSize / totalSize);
  VolumetricFeatures.push_back(neSize / totalSize);
  VolumetricFeatures.push_back(tuSize / totalSize);
  VolumetricFeatures.push_back(tuSize / neSize);
  VolumetricFeatures.push_back(edemaSize / (tuSize + neSize + edemaSize));
  
  return VolumetricFeatures;
}


VectorDouble ImagingSubtypePredictor::ScalingZeroToOne(const VectorDouble &data, const VariableLengthVectorType &minimum, const VariableLengthVectorType &maximum)
{
	VectorDouble scaledData;
	for (size_t i = 0; i < data.size(); i++)
		scaledData.push_back((data[i] - minimum[i]) / (maximum[i] - minimum[i]));
	return scaledData;

}

VectorDouble ImagingSubtypePredictor::CalculateEuclideanDistance(const VectorDouble &data, const VariableLengthVectorType &mean1, const VariableLengthVectorType &mean2, const VariableLengthVectorType &mean3)
{
	VectorDouble distances;
	float Dist1 = 0.0;
	float Dist2 = 0.0;
	float Dist3 = 0.0;
	for (size_t i = 0; i < data.size(); i++)
	{
		Dist1 += powf(data[i] - mean1[i], 2);
		Dist2 += powf(data[i] - mean2[i], 2);
		Dist3 += powf(data[i] - mean3[i], 2);
	}
	distances.push_back(sqrtf(Dist1));
	distances.push_back(sqrtf(Dist2));
	distances.push_back(sqrtf(Dist3));
	return distances;
}


