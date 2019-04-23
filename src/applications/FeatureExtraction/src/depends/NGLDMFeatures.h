/**
\file  NGLDMFeatures.h

\brief NGLDM feature calculation (Neighbouring Grey Level Dependence Matrix based features)

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

// ITK
#include "itkImage.h"
#include "itkHistogram.h"
#include "itkNumericTraits.h"
#include "itkVectorContainer.h"
#include "itkNeighborhoodIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMaskImageFilter.h"
#include "itkOffset.h"
#include "itkSize.h"
#include "itkImageRegionConstIterator.h"
//#include "itkEnhancedScalarImageToTextureFeaturesFilter.h"

// Eigen
#include "Eigen/Dense"

// STL
#include <sstream>
#include <map>
#include <string>
#include <numeric>
#include <vector>

// CaPTk
#include "FeatureBase.h"

//! Helper Class
struct NGLDMMatrixHolder
{
public:
	NGLDMMatrixHolder(double min, double max, int number, int depenence);

	int IntensityToIndex(double intensity);
	double IndexToMinIntensity(int index);
	double IndexToMeanIntensity(int index);
	double IndexToMaxIntensity(int index);

	double m_MinimumRange;
	double m_MaximumRange;
	double m_Stepsize;
	int m_NumberOfDependences;
	int m_NumberOfBins;
	Eigen::MatrixXd m_Matrix;

	int m_NeighbourhoodSize;
	unsigned long m_NumberOfNeighbourVoxels;
	unsigned long m_NumberOfDependenceNeighbourVoxels;
	unsigned long m_NumberOfNeighbourhoods;
	unsigned long m_NumberOfCompleteNeighbourhoods;
};

//! Helper Class
struct NGLDMMatrixFeatures
{
	NGLDMMatrixFeatures() :
		LowDependenceEmphasis(0),
		HighDependenceEmphasis(0),
		LowGreyLevelCountEmphasis(0),
		HighGreyLevelCountEmphasis(0),
		LowDependenceLowGreyLevelEmphasis(0),
		LowDependenceHighGreyLevelEmphasis(0),
		HighDependenceLowGreyLevelEmphasis(0),
		HighDependenceHighGreyLevelEmphasis(0),
		GreyLevelNonUniformity(0),
		GreyLevelNonUniformityNormalised(0),
		DependenceCountNonUniformity(0),
		DependenceCountNonUniformityNormalised(0),
		DependenceCountPercentage(0),
		GreyLevelVariance(0),
		DependenceCountVariance(0),
		DependenceCountEntropy(0),
		DependenceCountEnergy(0),
		MeanGreyLevelCount(0),
		MeanDependenceCount(0),
		ExpectedNeighbourhoodSize(0),
		AverageNeighbourhoodSize(0),
		AverageIncompleteNeighbourhoodSize(0),
		PercentageOfCompleteNeighbourhoods(0),
		PercentageOfDependenceNeighbours(0)
	{
	}

public:
	double LowDependenceEmphasis;
	double HighDependenceEmphasis;
	double LowGreyLevelCountEmphasis;
	double HighGreyLevelCountEmphasis;
	double LowDependenceLowGreyLevelEmphasis;
	double LowDependenceHighGreyLevelEmphasis;
	double HighDependenceLowGreyLevelEmphasis;
	double HighDependenceHighGreyLevelEmphasis;

	double GreyLevelNonUniformity;
	double GreyLevelNonUniformityNormalised;
	double DependenceCountNonUniformity;
	double DependenceCountNonUniformityNormalised;

	double DependenceCountPercentage;
	double GreyLevelVariance;
	double DependenceCountVariance;
	double DependenceCountEntropy;
	double DependenceCountEnergy;
	double MeanGreyLevelCount;
	double MeanDependenceCount;

	double ExpectedNeighbourhoodSize;
	double AverageNeighbourhoodSize;
	double AverageIncompleteNeighbourhoodSize;
	double PercentageOfCompleteNeighbourhoods;
	double PercentageOfDependenceNeighbours;
};


template< typename TImageType >
class NGLDMFeatures : public FeatureBase < TImageType >
{
public:
	//! Default constructor
	NGTDMFeatures() { };

	//! Default destructor
	~NGTDMFeatures() { };

	/**
	\brief Set the range/Chebyshev Distance (how far from the center index of interest do you want to calculate neighborhood level dependence); defaults to 1
	
	\IBSI calls this value Chebyshev Distance Delta.
	
	\param rangeValue integer value for the value you want for range
	**/
	void SetRange(int rangeValue)
	{
		m_range = rangeValue;
	}

	/**
	\brief Set the courseness/alpha parameter (How much difference between the grey levels of the neighbouring voxels do you consider to be dependent); defaults to 0 as typical choice for couraseness

	\IBSI calls this value alpha.

	\param coursenessValue integer value for the value you want for alpha
	**/
	void SetAlpha(int alphaValue)
	{
		m_alpha = alphaValue;
	}

	/**
	\brief Get the range/Chebyshev Distance (how far from the center index of interest do you want to calculate neighborhood level difference); returns int value
	**/
	int GetRange() {
		return m_range;
	}

	/**
	\brief Get the courseness/alpha parameter (How much difference between the grey levels of the neighbouring voxels do you consider to be dependent);  returns int value
	**/
	int GetAlpha() {
		return m_alpha;
	}

	/**
	\brief Get the number of bins for quantization of the input image
	**/
	int GetNumBins() {
		return m_bins;
	}

	/**
	\brief Set the number of bins for quantization of the input image; defaults to 10 unless overridden by user

	\param numBinValue integer value for the number of bins you want for image quantization
	**/
	void SetNumBins(int numBinValue)
	{
		m_bins = numBinValue;
	}

	/**
	\brief Set the minimum
	*/
	void SetMinimum(int minimumInput)
	{
		m_minimum = minimumInput;
	}

	/**
	\brief Set the maximum
	*/
	void SetMaximum(int maximumInput)
	{
		m_maximum = maximumInput;
	}

	/**
	\brief Update calculates features
	**/
	void Update()
	{
		if (!this->m_algorithmDone)
		{
			// do the computation that is needed here; update min max pixel values within ROI if not a
			if ((m_minimum == 0) || (m_maximum == 0))
			{
				//Masked out regions outside of the mask (ROI region, single label of multi label ROi).
				auto maskFilter = itk::MaskImageFilter< TImageType, TImageType, TImageType >::New();
				maskFilter->SetInput(this->m_inputImage); //full input image
				maskFilter->SetMaskImage(this->m_Mask);   //Assumed this is single label binary mask (already extracted from multi-label segmentation ROI)
				maskFilter->SetOutsideValue(0);
				maskFilter->Update();
				//Question: Is the mask being read in only one label of multi-label (e.g. edema of multi label segmentation?)

				//update m_minimum and m_maximum based on max and min pixel values of the ROI regions of the image
				auto minMaxComputer = itk::MinimumMaximumImageCalculator< TImageType >::New();
				minMaxComputer->SetImage(maskFilter->GetOutput());
				minMaxComputer->Compute();

				if (m_minimum == 0)
				{
					m_minimum = minMaxComputer->GetMinimum();
				}
				if (m_maximum == 0)
				{
					m_maximum = minMaxComputer->GetMaximum();
				}
			}

			std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - m_minimum = " << m_minimum << std::endl;
			std::cout << "\n[DEBUG] NGTDMFeatures.h - Update() - m_maximum = " << m_maximum << std::endl;

			typedef itk::NeighborhoodIterator<TImageType> ShapeIterType;
			typedef itk::NeighborhoodIterator<TImageType> ShapeMaskIterType;

			m_NumberOfCompleteNeighbourhoods = 0;
			m_NumberOfNeighbourhoods = 0;
			m_NumberOfNeighbourVoxels = 0;
			m_NumberOfDependenceNeighbourVoxels = 0;

			itk::Size<TImageType::ImageDimension> radius;
			radius.Fill(range);

			if ((direction > 1) && (direction - 2 < TImageType::ImageDimension))
			{
				radius[direction - 2] = 0;
			}

			ShapeIterType imageIter(radius, this->m_inputImage, this->m_inputImage->GetLargestPossibleRegion());
			ShapeMaskIterType maskIter(radius, this->m_Mask, this->m_Mask->GetLargestPossibleRegion());

			auto region = this->m_Mask->GetLargestPossibleRegion();

			auto center = imageIter.Size() / 2;
			auto iterSize = imageIter.Size();
			m_NeighbourhoodSize = iterSize - 1;
			while (!maskIter.IsAtEnd())
			{
				int sameValues = 0;
				bool completeNeighbourhood = true;

				int i = IntensityToIndex(imageIter.GetCenterPixel());

				if ((imageIter.GetCenterPixel() != imageIter.GetCenterPixel()) ||
					(maskIter.GetCenterPixel() < 1))
				{
					++imageIter;
					++maskIter;
					continue;
				}

				for (unsigned int position = 0; position < iterSize; ++position)
				{
					if (position == center)
					{
						continue;
					}
					if (!region.IsInside(maskIter.GetIndex(position)))
					{
						completeNeighbourhood = false;
						continue;
					}
					bool isInBounds;
					auto jIntensity = imageIter.GetPixel(position, isInBounds);
					auto jMask = maskIter.GetPixel(position, isInBounds);
					if (jMask < 1 || (jIntensity != jIntensity) || (!isInBounds))
					{
						completeNeighbourhood = false;
						continue;
					}

					int j = holder.IntensityToIndex(jIntensity);
					m_NumberOfNeighbourVoxels += 1;
					if (std::abs(i - j) <= alpha)
					{
						m_NumberOfDependenceNeighbourVoxels += 1;
						++sameValues;
					}
				}
				m_Matrix(i, sameValues) += 1;
				m_NumberOfNeighbourhoods += 1;
				if (completeNeighbourhood)
				{
					m_NumberOfCompleteNeighbourhoods += 1;
				}

				++imageIter;
				++maskIter;
			}

			this->m_algorithmDone = true;
		}
	};

	/**
	\brief return the map of feature names and feature values
	**/
	std::map< std::string, double > GetOutput()
	{
		if (!this->m_algorithmDone)
		{
			Update();
		}
		return this->m_features;
	};

	/**
	\brief return the Low Dependence Emphasis
	**/
	double GetLowDependenceEmphasis()
	{
		if (!this->m_algorithmDone)
		{
			Update();
		}
		return this->m_features["LowDependenceEmphasis"];
	}
	   
private:

	template<typename TPixel, unsigned int TImageType::Dimension>
	void
		CalculateNGLDMMatrix(itk::Image<TPixel, TImageType::Dimension>* itkImage,
			itk::Image<unsigned short, TImageType::Dimension>* mask,
			int alpha,
			int range,
			unsigned int direction,
			NGLDMMatrixHolder& holder)
	{
		typedef itk::Image<TPixel, TImageType::Dimension> ImageType;
		typedef itk::Image<unsigned short, TImageType::Dimension> MaskImageType;
		typedef itk::NeighborhoodIterator<ImageType> ShapeIterType;
		typedef itk::NeighborhoodIterator<MaskImageType> ShapeMaskIterType;

		holder.m_NumberOfCompleteNeighbourhoods = 0;
		holder.m_NumberOfNeighbourhoods = 0;
		holder.m_NumberOfNeighbourVoxels = 0;
		holder.m_NumberOfDependenceNeighbourVoxels = 0;

		itk::Size<TImageType::Dimension> radius;
		radius.Fill(range);

		if ((direction > 1) && (direction - 2 < TImageType::Dimension))
		{
			radius[direction - 2] = 0;
		}

		ShapeIterType imageIter(radius, itkImage, itkImage->GetLargestPossibleRegion());
		ShapeMaskIterType maskIter(radius, mask, mask->GetLargestPossibleRegion());

		auto region = mask->GetLargestPossibleRegion();

		auto center = imageIter.Size() / 2;
		auto iterSize = imageIter.Size();
		holder.m_NeighbourhoodSize = iterSize - 1;
		while (!maskIter.IsAtEnd())
		{
			int sameValues = 0;
			bool completeNeighbourhood = true;

			int i = holder.IntensityToIndex(imageIter.GetCenterPixel());

			if ((imageIter.GetCenterPixel() != imageIter.GetCenterPixel()) ||
				(maskIter.GetCenterPixel() < 1))
			{
				++imageIter;
				++maskIter;
				continue;
			}

			for (unsigned int position = 0; position < iterSize; ++position)
			{
				if (position == center)
				{
					continue;
				}
				if (!region.IsInside(maskIter.GetIndex(position)))
				{
					completeNeighbourhood = false;
					continue;
				}
				bool isInBounds;
				auto jIntensity = imageIter.GetPixel(position, isInBounds);
				auto jMask = maskIter.GetPixel(position, isInBounds);
				if (jMask < 1 || (jIntensity != jIntensity) || (!isInBounds))
				{
					completeNeighbourhood = false;
					continue;
				}

				int j = holder.IntensityToIndex(jIntensity);
				holder.m_NumberOfNeighbourVoxels += 1;
				if (std::abs(i - j) <= alpha)
				{
					holder.m_NumberOfDependenceNeighbourVoxels += 1;
					++sameValues;
				}
			}
			holder.m_Matrix(i, sameValues) += 1;
			holder.m_NumberOfNeighbourhoods += 1;
			if (completeNeighbourhood)
			{
				holder.m_NumberOfCompleteNeighbourhoods += 1;
			}

			++imageIter;
			++maskIter;
		}
	}


	void LocalCalculateFeatures(
		NGLDMMatrixHolder& holder,
		NGLDMMatrixFeatures& results
	)
	{
		auto sijMatrix = holder.m_Matrix;
		auto piMatrix = holder.m_Matrix;
		auto pjMatrix = holder.m_Matrix;

		double Ns = sijMatrix.sum();
		piMatrix.rowwise().normalize();
		pjMatrix.colwise().normalize();

		for (int i = 0; i < holder.m_NumberOfBins; ++i)
		{
			double sj = 0;
			for (int j = 0; j < holder.m_NumberOfDependences; ++j)
			{
				double iInt = i + 1;// holder.IndexToMeanIntensity(i);
				double sij = sijMatrix(i, j);
				double k = j + 1;
				double pij = sij / Ns;

				results.LowDependenceEmphasis += sij / k / k;
				results.HighDependenceEmphasis += sij * k * k;
				if (iInt != 0)
				{
					results.LowGreyLevelCountEmphasis += sij / iInt / iInt;
				}
				results.HighGreyLevelCountEmphasis += sij * iInt * iInt;
				if (iInt != 0)
				{
					results.LowDependenceLowGreyLevelEmphasis += sij / k / k / iInt / iInt;
				}
				results.LowDependenceHighGreyLevelEmphasis += sij * iInt * iInt / k / k;
				if (iInt != 0)
				{
					results.HighDependenceLowGreyLevelEmphasis += sij * k * k / iInt / iInt;
				}
				results.HighDependenceHighGreyLevelEmphasis += sij * k * k * iInt * iInt;


				results.MeanGreyLevelCount += iInt * pij;
				results.MeanDependenceCount += k * pij;
				if (pij > 0)
				{
					results.DependenceCountEntropy -= pij * std::log(pij) / std::log(2);
				}
				results.DependenceCountEnergy += pij * pij;
				sj += sij;
			}
			results.GreyLevelNonUniformity += sj * sj;
			results.GreyLevelNonUniformityNormalised += sj * sj;
		}

		for (int j = 0; j < holder.m_NumberOfDependences; ++j)
		{
			double si = 0;
			for (int i = 0; i < holder.m_NumberOfBins; ++i)
			{
				double sij = sijMatrix(i, j);
				si += sij;
			}
			results.DependenceCountNonUniformity += si * si;
			results.DependenceCountNonUniformityNormalised += si * si;
		}
		for (int i = 0; i < holder.m_NumberOfBins; ++i)
		{
			for (int j = 0; j < holder.m_NumberOfDependences; ++j)
			{
				double iInt = i + 1;// holder.IndexToMeanIntensity(i);
				double sij = sijMatrix(i, j);
				double k = j + 1;
				double pij = sij / Ns;

				results.GreyLevelVariance += (iInt - results.MeanGreyLevelCount) * (iInt - results.MeanGreyLevelCount) * pij;
				results.DependenceCountVariance += (k - results.MeanDependenceCount) * (k - results.MeanDependenceCount) * pij;
			}
		}
		results.LowDependenceEmphasis /= Ns;
		results.HighDependenceEmphasis /= Ns;
		results.LowGreyLevelCountEmphasis /= Ns;
		results.HighGreyLevelCountEmphasis /= Ns;
		results.LowDependenceLowGreyLevelEmphasis /= Ns;
		results.LowDependenceHighGreyLevelEmphasis /= Ns;
		results.HighDependenceLowGreyLevelEmphasis /= Ns;
		results.HighDependenceHighGreyLevelEmphasis /= Ns;

		results.GreyLevelNonUniformity /= Ns;
		results.GreyLevelNonUniformityNormalised /= (Ns * Ns);
		results.DependenceCountNonUniformity /= Ns;
		results.DependenceCountNonUniformityNormalised /= (Ns * Ns);

		results.DependenceCountPercentage = 1;

		results.ExpectedNeighbourhoodSize = holder.m_NeighbourhoodSize;
		results.AverageNeighbourhoodSize = holder.m_NumberOfNeighbourVoxels / (1.0 * holder.m_NumberOfNeighbourhoods);
		results.AverageIncompleteNeighbourhoodSize = (holder.m_NumberOfNeighbourVoxels - holder.m_NumberOfCompleteNeighbourhoods * holder.m_NeighbourhoodSize) / (1.0 * (holder.m_NumberOfNeighbourhoods - holder.m_NumberOfCompleteNeighbourhoods));
		results.PercentageOfCompleteNeighbourhoods = (1.0 * holder.m_NumberOfCompleteNeighbourhoods) / (1.0 * holder.m_NumberOfNeighbourhoods);
		results.PercentageOfDependenceNeighbours = holder.m_NumberOfDependenceNeighbourVoxels / (1.0 * holder.m_NumberOfNeighbourVoxels);
	
		std::cout << "[DEBUG] LowDependenceEmphasis = " << results.LowDependenceEmphasis << std::endl;
		this->m_features["LowDependenceEmphasis"] = results.LowDependenceEmphasis;

		std::cout << "[DEBUG] HighDependenceEmphasis = " << results.HighDependenceEmphasis << std::endl;
		this->m_features["HighDependenceEmphasis"] = results.HighDependenceEmphasis;
	}

	//template<typename TPixel, unsigned int TImageType::Dimension>
	//void
	//	CalculateCoocurenceFeatures(itk::Image<TPixel, TImageType::Dimension>* itkImage, TImageType::Pointer mask)
	//{
	//	typedef itk::Image<unsigned short, TImageType::Dimension> MaskType;

	//	double rangeMin = config.MinimumIntensity;
	//	double rangeMax = config.MaximumIntensity;
	//	int numberOfBins = config.Bins;

	//	typename MaskType::Pointer maskImage = MaskType::New();
	//	mitk::CastToItkImage(mask, maskImage);

	//	std::vector<NGLDMMatrixFeatures> resultVector;
	//	int numberofDependency = 37;
	//	if (TImageType::Dimension == 2)
	//		numberofDependency = 37;

	//	NGLDMMatrixHolder holderOverall(rangeMin, rangeMax, numberOfBins, numberofDependency);
	//	NGLDMMatrixFeatures overallFeature;
	//	CalculateNGLDMMatrix<TPixel, TImageType::Dimension>(itkImage, maskImage, config.alpha, config.range, config.direction, holderOverall);
	//	LocalCalculateFeatures(holderOverall, overallFeature);

	//	//MatrixFeaturesTo(overallFeature, config.FeatureEncoding, featureList);
	//}


	//static
	//	void MatrixFeaturesTo(NGLDMMatrixFeatures features,
	//		std::string prefix,
	//		GIFNeighbouringGreyLevelDependenceFeature::FeatureListType & featureList)
	//{
	//	featureList.push_back(std::make_pair(prefix + "Low Dependence Emphasis", features.LowDependenceEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "High Dependence Emphasis", features.HighDependenceEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "Low Grey Level Count Emphasis", features.LowGreyLevelCountEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "High Grey Level Count Emphasis", features.HighGreyLevelCountEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "Low Dependence Low Grey Level Emphasis", features.LowDependenceLowGreyLevelEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "Low Dependence High Grey Level Emphasis", features.LowDependenceHighGreyLevelEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "High Dependence Low Grey Level Emphasis", features.HighDependenceLowGreyLevelEmphasis));
	//	featureList.push_back(std::make_pair(prefix + "High Dependence High Grey Level Emphasis", features.HighDependenceHighGreyLevelEmphasis));

	//	featureList.push_back(std::make_pair(prefix + "Grey Level Non-Uniformity", features.GreyLevelNonUniformity));
	//	featureList.push_back(std::make_pair(prefix + "Grey Level Non-Uniformity Normalised", features.GreyLevelNonUniformityNormalised));
	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Non-Uniformity", features.DependenceCountNonUniformity));
	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Non-Uniformity Normalised", features.DependenceCountNonUniformityNormalised));

	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Percentage", features.DependenceCountPercentage));
	//	featureList.push_back(std::make_pair(prefix + "Grey Level Mean", features.MeanGreyLevelCount));
	//	featureList.push_back(std::make_pair(prefix + "Grey Level Variance", features.GreyLevelVariance));
	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Mean", features.MeanDependenceCount));
	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Variance", features.DependenceCountVariance));
	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Entropy", features.DependenceCountEntropy));
	//	featureList.push_back(std::make_pair(prefix + "Dependence Count Energy", features.DependenceCountEnergy));

	//	featureList.push_back(std::make_pair(prefix + "Expected Neighbourhood Size", features.ExpectedNeighbourhoodSize));
	//	featureList.push_back(std::make_pair(prefix + "Average Neighbourhood Size", features.AverageNeighbourhoodSize));
	//	featureList.push_back(std::make_pair(prefix + "Average Incomplete Neighbourhood Size", features.AverageIncompleteNeighbourhoodSize));
	//	featureList.push_back(std::make_pair(prefix + "Percentage of complete Neighbourhoods", features.PercentageOfCompleteNeighbourhoods));
	//	featureList.push_back(std::make_pair(prefix + "Percentage of Dependence Neighbour Voxels", features.PercentageOfDependenceNeighbours));

	//}

	//simple functions
	int NGLDMMatrixHolder::IntensityToIndex(double intensity)
	{
		return std::floor((intensity - m_MinimumRange) / m_Stepsize);
	}

	double NGLDMMatrixHolder::IndexToMinIntensity(int index)
	{
		return m_MinimumRange + index * m_Stepsize;
	}

	double NGLDMMatrixHolder::IndexToMeanIntensity(int index)
	{
		return m_MinimumRange + (index + 0.5) * m_Stepsize;
	}

	double NGLDMMatrixHolder::IndexToMaxIntensity(int index)
	{
		return m_MinimumRange + (index + 1) * m_Stepsize;
	}

	//variables
	unsigned int m_range = 1;
	unsigned int m_bins = 10;
	typename TImageType::PixelType m_minimum = 0;
	typename TImageType::PixelType m_maximum = 0;

	typename TImageType::SizeType m_radius;
		
	double m_MinimumRange;
	double m_MaximumRange;
	double m_Stepsize;
	int m_NumberOfDependences;
	Eigen::MatrixXd m_Matrix;

	int m_NeighbourhoodSize;
	unsigned long m_NumberOfNeighbourVoxels;
	unsigned long m_NumberOfDependenceNeighbourVoxels;
	unsigned long m_NumberOfNeighbourhoods;
	unsigned long m_NumberOfCompleteNeighbourhoods;
};