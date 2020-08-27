/**
\file  NGLDMFeatures.h

\brief NGLDM feature calculation (Neighbouring Grey Level Dependence Matrix based features)

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html

*/
#pragma once

// ITK
#include "itkImage.h"
//#include "itkHistogram.h"
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
#include "TextureFeatureBase.h"

namespace mitk
{
	//! Helper Class
	struct NGLDMMatrixHolder
	{
	public:
		NGLDMMatrixHolder(double min, double max, int number, int depenence);

		int IntensityToIndex(double intensity);
		double IndexToMinIntensity(int index);
		double IndexToMeanIntensity(int index);
		double IndexToMaxIntensity(int index);
		void Print();

		double m_minimumRange;
		double m_maximumRange;
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
}

mitk::NGLDMMatrixHolder::NGLDMMatrixHolder(double min, double max, int number, int depenence) :
	m_minimumRange(min),
	m_maximumRange(max),
	m_Stepsize(0),
	m_NumberOfDependences(depenence),
	m_NumberOfBins(number),
	m_NeighbourhoodSize(1),
	m_NumberOfNeighbourVoxels(0),
	m_NumberOfDependenceNeighbourVoxels(0),
	m_NumberOfNeighbourhoods(0),
	m_NumberOfCompleteNeighbourhoods(0)
{
	m_Matrix.resize(number, depenence);
	//std::cout << "[DEBUG] NGLDMFeatures.h::mitk::NGLDMMatrixHolder::NGLDMMatrixHolder:: number = " << number << " | depenence = " << depenence << std::endl;

	m_Matrix.fill(0);
	//std::cout << "[DEBUG] NGLDMFeatures.h::mitk::NGLDMMatrixHolder::NGLDMMatrixHolder:: m_Matrix = " << m_Matrix << std::endl;

	m_Stepsize = (max - min) / (number);
}


void mitk::NGLDMMatrixHolder::Print()
{
	//std::cout << "[DEBUG] NGLDMFeatures.h::mitk::NGLDMMatrixHolder::NGLDMMatrixHolder::Print:: m_Matrix = \n\n" << m_Matrix << std::endl;
}

int mitk::NGLDMMatrixHolder::IntensityToIndex(double intensity)
{
	return std::floor((intensity - m_minimumRange) / m_Stepsize);
}

double mitk::NGLDMMatrixHolder::IndexToMinIntensity(int index)
{
	return m_minimumRange + index * m_Stepsize;
}
double mitk::NGLDMMatrixHolder::IndexToMeanIntensity(int index)
{
	return m_minimumRange + (index + 0.5) * m_Stepsize;
}
double mitk::NGLDMMatrixHolder::IndexToMaxIntensity(int index)
{
	return m_minimumRange + (index + 1) * m_Stepsize;
}

template< typename TImageType >
class NGLDMFeatures : 
  public FeatureBase < TImageType >, public TextureFeatureBase< TImageType >
{
public:
	//! Default constructor
	NGLDMFeatures() { };

	//! Default destructor
	~NGLDMFeatures() { };

  void SetDistanceMax(const double maxDistance)
  {
    m_maxDistance = maxDistance;
  }

	/**
	\brief Set the courseness/alpha parameter (How much difference between the grey levels of the neighbouring voxels do you consider to be dependent); defaults to 0 as typical choice for couraseness

	\IBSI calls this value alpha.

	\param coursenessValue integer value for the value you want for alpha
	**/
	// void SetAlpha(int alphaValue)
	// {
	// 	m_alpha = alphaValue;
	// }

	/**
	\brief Update calculates features
	**/
	void Update()
	{
		if (!this->m_algorithmDone)
		{
			// do the computation that is needed here; update min max pixel values within ROI if not a
			if ((this->m_minimum == 0) || (this->m_maximum == 0))
			{
				//Masked out regions outside of the mask (ROI region, single label of multi label ROi).
				auto maskFilter = itk::MaskImageFilter< TImageType, TImageType, TImageType >::New();
				maskFilter->SetInput(this->m_inputImage); //full input image
				maskFilter->SetMaskImage(this->m_Mask);   //Assumed this is single label binary mask (already extracted from multi-label segmentation ROI)
				maskFilter->SetOutsideValue(0);
				maskFilter->Update();
				//Question: Is the mask being read in only one label of multi-label (e.g. edema of multi label segmentation?)

				//update this->m_minimum and this->m_maximum based on max and min pixel values of the ROI regions of the image
				auto minMaxComputer = itk::MinimumMaximumImageCalculator< TImageType >::New();
				minMaxComputer->SetImage(maskFilter->GetOutput());
				minMaxComputer->Compute();

				if (this->m_minimum == 0)
				{
					this->m_minimum = minMaxComputer->GetMinimum();
				}
				if (this->m_maximum == 0)
				{
					this->m_maximum = minMaxComputer->GetMaximum();
				}
			}

			//std::cout << "\n[DEBUG] NGLDMFeatures.h - Update() - this->m_minimum = " << this->m_minimum << std::endl;
			//std::cout << "\n[DEBUG] NGLDMFeatures.h - Update() - this->m_maximum = " << this->m_maximum << std::endl;

			//TBD - Hard Coded for Phantom for now
			int alpha = 0;
			unsigned int direction = 26;
			int max_number_of_dependence = 21;
			int range = 1;
			//TBD - Hard Coded for Phantom for now


			//NGLDMMatrixHolder holderOverall(rangeMin, rangeMax, numberOfBins, numberofDependency);
			//std::cout << "[DEBUG] mitk::NGLDMMatrixHolder::holderOverall(" << this->m_minimum << ", " << this->m_maximum << ", " << this->m_Bins << ", " << max_number_of_dependence << ")" << std::endl;
			mitk::NGLDMMatrixHolder holderOverall(this->m_minimum, this->m_maximum, this->m_Bins, max_number_of_dependence);
			mitk::NGLDMMatrixFeatures overallFeature;
			CalculateNGLDMMatrix(this->m_inputImage, this->m_Mask, alpha, range, direction, holderOverall);
			//std::cout << "[DEBUG] NGLDMFeatures.h||print " << std::endl;
			holderOverall.Print();

			LocalCalculateFeatures(holderOverall, overallFeature);

		}
	};

private:

	void CalculateNGLDMMatrix(
		typename TImageType::Pointer itkImage,
		typename TImageType::Pointer mask,
		int alpha,
		int range,
		unsigned int direction,
		mitk::NGLDMMatrixHolder& holder
	)
	{
		// typedef itk::Image< TImageType > ImageType;
		// typedef itk::Image< TImageType > MaskImageType;
		typedef itk::NeighborhoodIterator<TImageType> ShapeIterType;
		typedef itk::NeighborhoodIterator<TImageType> ShapeMaskIterType;

		holder.m_NumberOfCompleteNeighbourhoods = 0;
		holder.m_NumberOfNeighbourhoods = 0;
		holder.m_NumberOfNeighbourVoxels = 0;
		holder.m_NumberOfDependenceNeighbourVoxels = 0;

		itk::Size<TImageType::ImageDimension> radius;
		radius.Fill(range);

		if ((direction > 1) && (direction - 2 < TImageType::ImageDimension))
		{
			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::direction = " << direction << std::endl;
			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::TImageType::ImageDimension = " << TImageType::ImageDimension << std::endl;
			radius[direction - 2] = 0;
		}

		//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::radius = " << radius << std::endl;

		ShapeIterType imageIter(radius, itkImage, itkImage->GetLargestPossibleRegion());
		ShapeMaskIterType maskIter(radius, mask, mask->GetLargestPossibleRegion());

		auto region = mask->GetLargestPossibleRegion();

		auto center = imageIter.Size() / 2;
		//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::center = " << center << std::endl; //phantom 13

		auto iterSize = imageIter.Size();
		//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::iterSize = " << iterSize << std::endl; //phantom 27

		holder.m_NeighbourhoodSize = iterSize - 1; //phantom 26
		
		int iteration_count = 0;
		while (!maskIter.IsAtEnd())
		{
			int sameValues = 0;
			bool completeNeighbourhood = true;

			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::imageIter.GetCenterPixel() = " << imageIter.GetCenterPixel() << std::endl;

			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::imageIter.GetCenterPointer() = " << imageIter.GetCenterPointer() << std::endl;

			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::iteration_count = " << iteration_count << "\n\n"<< std::endl;

			int i = holder.IntensityToIndex(imageIter.GetCenterPixel());
			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << std::endl;
			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i intensity = " << imageIter.GetCenterPixel() << std::endl;


			//if ((imageIter.GetCenterPixel() != imageIter.GetCenterPixel()) || (maskIter.GetCenterPixel() < 1))
			if ((imageIter.GetCenterPixel() != maskIter.GetCenterPixel()) || (maskIter.GetCenterPixel() < 1))
			{
				//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | outside ROI." << std::endl;

				++imageIter;
				++maskIter;
				continue;
			}

			for (unsigned int position = 0; position < iterSize; ++position)
			{
				//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | position = " << position << "/" << iterSize << std::endl;

				if (position == center)
				{
					//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | position = " << position << "/" << iterSize << " == center = " << center << " | Skip" << std::endl;

					continue;
				}
				if (!region.IsInside(maskIter.GetIndex(position)))
				{
					completeNeighbourhood = false;
					continue;
				}
				else {
					//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | position = " << position << "/" << iterSize << " | completeNeighbourhood == " << completeNeighbourhood << std::endl;
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

				//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::j (Index of CenterVoxel) = " << j << std::endl;
				//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::j intensity = " << jIntensity << std::endl;

				holder.m_NumberOfNeighbourVoxels += 1;
				//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | m_NumberOfNeighbourVoxels = " << holder.m_NumberOfNeighbourVoxels << std::endl;

				if (std::abs(i - j) <= alpha)
				{
					holder.m_NumberOfDependenceNeighbourVoxels += 1;
					++sameValues;
				}
			}

			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | m_NumberOfNeighbourVoxels = " << holder.m_NumberOfNeighbourVoxels << std::endl;

			holder.m_Matrix(i, sameValues) += 1;
			holder.m_NumberOfNeighbourhoods += 1;

			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << "| m_Matrix = " << holder.m_Matrix << std::endl;

			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | m_Matrix(" << i << ", " << sameValues << ") = " << holder.m_Matrix(i, sameValues) << std::endl;
			//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | m_NumberOfNeighbourhoods = " << holder.m_NumberOfNeighbourhoods << std::endl;

			if (completeNeighbourhood)
			{
				holder.m_NumberOfCompleteNeighbourhoods += 1;
				//std::cout << "[DEBUG] NGLDMFeatures.h::CalculateNGLDMMatrix::maskInteration::i (Index of CenterVoxel) = " << i << " | completeNeighbourhood = " << completeNeighbourhood << std::endl;
			}

			++imageIter;
			++maskIter;

			//TBD
			++iteration_count;
		}

		//std::cout << "[DEBUG] m_NumberOfCompleteNeighbourhoods = " << holder.m_NumberOfCompleteNeighbourhoods << ::std::endl;
		//std::cout << "[DEBUG] m_NumberOfNeighbourhoods = " << holder.m_NumberOfNeighbourhoods << ::std::endl;
		//std::cout << "[DEBUG] m_NumberOfNeighbourVoxels = " << holder.m_NumberOfNeighbourVoxels << ::std::endl;
		//std::cout << "[DEBUG] m_NumberOfDependenceNeighbourVoxels = " << holder.m_NumberOfDependenceNeighbourVoxels << ::std::endl;
	}


	void LocalCalculateFeatures(
		mitk::NGLDMMatrixHolder& holder,
		mitk::NGLDMMatrixFeatures& results
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

		//std::cout << "[DEBUG] LowDependenceEmphasis = " << results.LowDependenceEmphasis << std::endl;
		this->m_features["LowDependenceEmphasis"] = results.LowDependenceEmphasis;

		//std::cout << "[DEBUG] HighDependenceEmphasis = " << results.HighDependenceEmphasis << std::endl;
		this->m_features["HighDependenceEmphasis"] = results.HighDependenceEmphasis;

		//std::cout << "[DEBUG] LowGreyLevelCountEmphasis = " << results.LowGreyLevelCountEmphasis << std::endl;
		this->m_features["LowGreyLevelCountEmphasis"] = results.LowGreyLevelCountEmphasis;

		//std::cout << "[DEBUG] HighGreyLevelCountEmphasis = " << results.HighGreyLevelCountEmphasis << std::endl;
		this->m_features["HighGreyLevelCountEmphasis"] = results.HighGreyLevelCountEmphasis;

		//std::cout << "[DEBUG] LowDependenceLowGreyLevelEmphasis = " << results.LowDependenceLowGreyLevelEmphasis << std::endl;
		this->m_features["LowDependenceLowGreyLevelEmphasis"] = results.LowDependenceLowGreyLevelEmphasis;

		//std::cout << "[DEBUG] LowDependenceHighGreyLevelEmphasis = " << results.LowDependenceHighGreyLevelEmphasis << std::endl;
		this->m_features["LowDependenceHighGreyLevelEmphasis"] = results.LowDependenceHighGreyLevelEmphasis;

		//std::cout << "[DEBUG] HighDependenceLowGreyLevelEmphasis = " << results.HighDependenceLowGreyLevelEmphasis << std::endl;
		this->m_features["HighDependenceLowGreyLevelEmphasis"] = results.HighDependenceLowGreyLevelEmphasis;

		//std::cout << "[DEBUG] HighDependenceHighGreyLevelEmphasis = " << results.HighDependenceHighGreyLevelEmphasis << std::endl;
		this->m_features["HighDependenceHighGreyLevelEmphasis"] = results.HighDependenceHighGreyLevelEmphasis;

		//std::cout << "[DEBUG] GreyLevelNonUniformity = " << results.GreyLevelNonUniformity << std::endl;
		this->m_features["GreyLevelNonUniformity"] = results.GreyLevelNonUniformity;

		//std::cout << "[DEBUG] GreyLevelNonUniformityNormalised = " << results.GreyLevelNonUniformityNormalised << std::endl;
		this->m_features["GreyLevelNonUniformityNormalised"] = results.GreyLevelNonUniformityNormalised;

		//std::cout << "[DEBUG] DependenceCountNonUniformity = " << results.DependenceCountNonUniformity << std::endl;
		this->m_features["DependenceCountNonUniformity"] = results.DependenceCountNonUniformity;

		//std::cout << "[DEBUG] DependenceCountNonUniformityNormalised = " << results.DependenceCountNonUniformityNormalised << std::endl;
		this->m_features["DependenceCountNonUniformityNormalised"] = results.DependenceCountNonUniformityNormalised;

		//std::cout << "[DEBUG] DependenceCountPercentage = " << results.DependenceCountPercentage << std::endl;
		this->m_features["DependenceCountPercentage"] = results.DependenceCountPercentage;

		//std::cout << "[DEBUG] GreyLevelVariance = " << results.GreyLevelVariance << std::endl;
		this->m_features["GreyLevelVariance"] = results.GreyLevelVariance;

		//std::cout << "[DEBUG] DependenceCountVariance = " << results.DependenceCountVariance << std::endl;
		this->m_features["DependenceCountVariance"] = results.DependenceCountVariance;
		
		//std::cout << "[DEBUG] DependenceCountEntropy = " << results.DependenceCountEntropy << std::endl;
		this->m_features["DependenceCountEntropy"] = results.DependenceCountEntropy;

		//std::cout << "[DEBUG] DependenceCountEnergy = " << results.DependenceCountEnergy << std::endl;
		this->m_features["DependenceCountEnergy"] = results.DependenceCountEnergy;

		//std::cout << "[DEBUG] MeanGreyLevelCount = " << results.MeanGreyLevelCount << std::endl;
		this->m_features["MeanGreyLevelCount"] = results.MeanGreyLevelCount;

		//std::cout << "[DEBUG] MeanDependenceCount = " << results.MeanDependenceCount << std::endl;
		this->m_features["MeanDependenceCount"] = results.MeanDependenceCount;


		//TBD - Extra feature values not listed in IBSI
		//double ExpectedNeighbourhoodSize;
		//double AverageNeighbourhoodSize;
		//double AverageIncompleteNeighbourhoodSize;
		//double PercentageOfCompleteNeighbourhoods;
		//double PercentageOfDependenceNeighbours;

	}

	//double m_minimumRange;
	//double m_maximumRange;
	//double m_Stepsize;
	//int m_NumberOfDependences;
	//Eigen::MatrixXd m_Matrix;

	//int m_NeighbourhoodSize;
	//unsigned long m_NumberOfNeighbourVoxels;
	//unsigned long m_NumberOfDependenceNeighbourVoxels;
	//unsigned long m_NumberOfNeighbourhoods;
	//unsigned long m_NumberOfCompleteNeighbourhoods;

	private:
    double m_maxDistance = -1;

    //variables
    unsigned int m_range = 1;
    typename TImageType::SizeType m_radius;

};