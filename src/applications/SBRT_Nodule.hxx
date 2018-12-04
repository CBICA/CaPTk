#pragma once

#include <sstream>
#include "SBRT_Nodule.h"

//using namespace cv;

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::SBRT_Nodule()
{
	for (int i=0; i<TInputImageNum; i++)
		m_modailityWt.push_back(0.5);

    m_lf_lab = 2;
	m_disThreshold = 5.0;
	m_sigma = 1.0;
	m_suvRatio = 3.0;
	m_iterNum = 200;
	m_smoothIterNum = 3;
	m_smoothR = 1;
	m_minNumFgSeed = 25;

	m_nei_r = 1;
	m_alpha = 0.9;

	m_outputBaseName = "seg";

	m_logger = NULL;
	m_logger_or_not = false;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::~SBRT_Nodule()
{

}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
int SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::Initialization( std::vector<std::string> &InputImageName, std::string maskName )
{
	if (m_logger_or_not)
		m_logger->Write("Read images ...");
	else
		std::cout << "Read images ... " << std::endl;	
	
	typename ScalarImageReaderType::Pointer InputReader = ScalarImageReaderType::New();
	
	typename ScalarMedianFilterType::Pointer medianFilter = ScalarMedianFilterType::New();
	typename ScalarMedianFilterType::InputSizeType mfRadius;
	mfRadius.Fill(m_smoothR);
	medianFilter->SetRadius(mfRadius);

	typename ScalarImageReaderType::Pointer MaskReader = ScalarImageReaderType::New();
	MaskReader->SetFileName( maskName );
	MaskReader->Update();

	for ( unsigned int i=0;i<TInputImageNum;i++ )
	{
		InputReader->SetFileName( InputImageName[i] );

		try
		{
			InputReader->Update();
		}
		catch ( itk::ExceptionObject &err )
		{
			std::cerr << " ExceptionObject caught !" << std::endl;
			std::cerr << " " << err << std::endl;
			return EXIT_FAILURE;
		}
		
		medianFilter->SetInput(InputReader->GetOutput());
		medianFilter->Update();

		ScalarImageTypePointer tmpImage = ScalarImageType::New();
		tmpImage->SetOrigin( InputReader->GetOutput()->GetOrigin() );
		tmpImage->SetSpacing( InputReader->GetOutput()->GetSpacing() );
		tmpImage->SetRegions( InputReader->GetOutput()->GetRequestedRegion() );
		tmpImage->SetDirection( InputReader->GetOutput()->GetDirection() );
		tmpImage->Allocate(); 
		tmpImage->FillBuffer( -1 );

		ScalarIteratorType MFiter( medianFilter->GetOutput(), medianFilter->GetOutput()->GetRequestedRegion() );
		ScalarIteratorType tmpiter( tmpImage, tmpImage->GetRequestedRegion() );
		MFiter.GoToBegin();
		tmpiter.GoToBegin();
		while ( !MFiter.IsAtEnd() )
		{
			tmpiter.Set( MFiter.Get() );

			++MFiter;
			++tmpiter;
		}

		m_MfImage.push_back(tmpImage);
	}

	m_InputImageVector = VectorImageType::New();
	m_InputImageVector->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_InputImageVector->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_InputImageVector->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_InputImageVector->SetDirection( m_MfImage[0]->GetDirection() );
	m_InputImageVector->Allocate();
	VectorType temp(0.0); 
	m_InputImageVector->FillBuffer( temp );

	VectorIteratorType Viter( m_InputImageVector,m_InputImageVector->GetRequestedRegion() );
	ScalarIteratorType Miter( MaskReader->GetOutput(), MaskReader->GetOutput()->GetRequestedRegion() );
	for ( unsigned int i=0;i<TInputImageNum;i++ )
	{
		ScalarIteratorType Siter( m_MfImage[i],m_MfImage[i]->GetRequestedRegion() );
		Siter.GoToBegin();
		Viter.GoToBegin();
		Miter.GoToBegin();
		while( !Viter.IsAtEnd() )
		{
			if ( Miter.Get()==m_lf_lab )
			{
				temp = Viter.Get();
				temp[i] = Siter.Get();
				Viter.Set(temp);
			}
			else
				Siter.Set(0.0);
			
			++Viter;
			++Siter;
			++Miter;
		}
	}

	// calculate image gradient
	m_GraMagnitude = ScalarImageType::New();
	m_GraMagnitude->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_GraMagnitude->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_GraMagnitude->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_GraMagnitude->SetDirection( m_MfImage[0]->GetDirection() );
	m_GraMagnitude->Allocate(); 
	m_GraMagnitude->FillBuffer( -1 );

	if ( 1==TInputImageNum )
	{
		typename ScalarGraMagFilterType::Pointer gmFilter = ScalarGraMagFilterType::New();
		gmFilter->SetInput(m_MfImage[0]);

		ScalarIteratorType Giter( gmFilter->GetOutput(), gmFilter->GetOutput()->GetRequestedRegion() );
		ScalarIteratorType MGiter( m_GraMagnitude, m_GraMagnitude->GetRequestedRegion() );
		Giter.GoToBegin();
		MGiter.GoToBegin();
		while ( !Giter.IsAtEnd() )
		{
			MGiter.Set( Giter.Get() );

			++Giter;
			++MGiter;
		}		
	}
	else if ( 1<TInputImageNum )
	{
		typename VectorGraMagFilterType::Pointer gmFilter = VectorGraMagFilterType::New();
		gmFilter->SetInput(m_InputImageVector);
		//parameters for the magnitude computation
		gmFilter->SetUsePrincipleComponentsOff();
		gmFilter->SetUseImageSpacingOn();

		ScalarIteratorType Giter( gmFilter->GetOutput(), gmFilter->GetOutput()->GetRequestedRegion() );
		ScalarIteratorType MGiter( m_GraMagnitude, m_GraMagnitude->GetRequestedRegion() );
		Giter.GoToBegin();
		MGiter.GoToBegin();
		while ( !Giter.IsAtEnd() )
		{
			MGiter.Set( Giter.Get() );

			++Giter;
			++MGiter;
		}
	}

	// initialize output image
	m_OutputImage = ScalarImageType::New();
	m_OutputImage->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_OutputImage->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_OutputImage->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_OutputImage->SetDirection( m_MfImage[0]->GetDirection() );
	m_OutputImage->Allocate(); 
	m_OutputImage->FillBuffer( -1 );
	// initialize valid voxel ID image
	m_voxelID = ScalarImageType::New();
	m_voxelID->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_voxelID->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_voxelID->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_voxelID->SetDirection( m_MfImage[0]->GetDirection() );
	m_voxelID->Allocate(); 
	m_voxelID->FillBuffer( 0 );

	ScalarIteratorType Iiter( m_voxelID, m_voxelID->GetRequestedRegion() );
	Iiter.GoToBegin();
	Miter.GoToBegin();
	m_validVxNum = 0;
	while ( !Miter.IsAtEnd() )
	{		
		if ( Miter.Get()==m_lf_lab )
		{
			m_validVxNum++;
			Iiter.Set( m_validVxNum );
		}
        else if ( Miter.Get()>0 )
        {
            Iiter.Set( -1 );
        }

		++Iiter;
		++Miter;
	}

	//m_MaskImage = MaskReader->GetOutput();
	//m_MaskImage->DisconnectPipeline();

	for( unsigned int j=0;j<TDimension;j++ )
	{
		m_ImageVoxelSize.push_back( m_InputImageVector->GetSpacing()[j] );
	}

	m_ImageSize = m_InputImageVector->GetRequestedRegion().GetSize();

	//if (m_logger_or_not)
	//	m_logger->Write("Read images successfully ...");
	//else
	//	cout << "Read images successfully " << endl;
    
    return 0;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::SetParameters( int lfLab, std::string obasename, float suvR, int mnfs, std::vector<float> &mWt, float dThres, float wSigma, int iter, int siter, int sr )
{
	for (int i=0; i<TInputImageNum; i++)
		m_modailityWt[i] = mWt[i];

    m_lf_lab = lfLab;
	m_disThreshold = dThres;
	m_sigma = wSigma;
	m_suvRatio = suvR;
	m_iterNum = iter;
	m_smoothIterNum = siter;
	m_smoothR = sr;
	m_minNumFgSeed = mnfs;

	m_outputBaseName = obasename;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::SetLogger(std::string logName)
{
	m_logger = new cbica::Logging;
	m_logger->UseNewFile(logName);
	m_logger_or_not = true;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::ReadSeedImage( std::string seedName )
{
	if (m_logger_or_not)
		m_logger->Write("Read seeds ...");
	else
		std::cout << "Read seeds..." << std::endl;
    
	typename ScalarImageReaderType::Pointer SeedReader = ScalarImageReaderType::New();
    SeedReader->SetFileName(seedName);
    SeedReader->Update();
    
    m_SeedImage = SeedReader->GetOutput();
    m_SeedImage->DisconnectPipeline();
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::GenerateSeeds()
{
  if (m_logger_or_not)
    m_logger->Write("Generate seeds ...");
  else
    std::cout << "Generate seeds..." << std::endl;
	
	m_SeedImage = ScalarImageType::New();
	m_SeedImage->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_SeedImage->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_SeedImage->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_SeedImage->SetDirection( m_MfImage[0]->GetDirection() );
	m_SeedImage->Allocate(); 
	m_SeedImage->FillBuffer( 0.0 );

	typename ScalarMinMaxImageCalculatorType::Pointer imageMinMaxCalc = ScalarMinMaxImageCalculatorType::New();
	imageMinMaxCalc->SetImage(m_MfImage[0]);
	imageMinMaxCalc->Compute();

	float maxSuv = imageMinMaxCalc->GetMaximum();
	float ratioSuv = maxSuv / m_suvRatio;

	//for test
	//cout << maxSuv << endl;

	typename ScalarThresholdImageFilterType::Pointer thresholdImage = ScalarThresholdImageFilterType::New();
	thresholdImage->SetInput(m_MfImage[0]);
	thresholdImage->ThresholdOutside(ratioSuv, maxSuv);
	thresholdImage->SetOutsideValue(0);
	thresholdImage->Update();

	m_labelingImage = ScalarLabelImageType::New();
	m_labelingImage->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_labelingImage->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_labelingImage->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_labelingImage->SetDirection( m_MfImage[0]->GetDirection() );
	m_labelingImage->Allocate(); 
	m_labelingImage->FillBuffer( 0 );

	ScalarLabelImageIteratorType lliter(m_labelingImage, m_labelingImage->GetRequestedRegion());
	ScalarIteratorType tiiter(thresholdImage->GetOutput(), thresholdImage->GetOutput()->GetRequestedRegion());
	lliter.GoToBegin();
	tiiter.GoToBegin();
	while( !tiiter.IsAtEnd() )
	{
		if (tiiter.Get()!=0)
			lliter.Set(1);

		++lliter;
		++tiiter;
	}

	typename ScalarConnectedComponentImageFilterType::Pointer connectedComp = ScalarConnectedComponentImageFilterType::New();
	connectedComp->SetInput(m_labelingImage);
	connectedComp->Update();

	int numCc = connectedComp->GetObjectCount();
	std::vector<float> tmpSuv;
	std::vector<std::vector<float> > ccSuvSet(numCc, tmpSuv);

	//for test
	//cout << numCc << endl;
	
	ScalarLabelImageIteratorType nciter( connectedComp->GetOutput(), connectedComp->GetOutput()->GetRequestedRegion() );
	ScalarIteratorType iiter( m_MfImage[0], m_MfImage[0]->GetRequestedRegion() );
	ScalarIteratorType sditer( m_SeedImage, m_SeedImage->GetRequestedRegion() );

	nciter.GoToBegin();
	iiter.GoToBegin();
	sditer.GoToBegin();
	while( !nciter.IsAtEnd() )
	{
		if (nciter.Get()>0)
		{
			ccSuvSet[nciter.Get()-1].push_back( iiter.Get() );
			sditer.Set(1);
		}

		++nciter;
		++iiter;
		++sditer;
	}
	// remove very small regions
	int smallNum = 5;

	nciter.GoToBegin();
	sditer.GoToBegin();
	while( !nciter.IsAtEnd() )
	{
		if (nciter.Get()>0)
		{
			if (ccSuvSet[nciter.Get()-1].size()<=smallNum)
				sditer.Set(0);
		}

		++nciter;
		++sditer;
	}

	// get background seeds, labeling the as 1
	StructureElementType structureElement;
	structureElement.SetRadius(5);
	structureElement.CreateStructuringElement();

	typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput(m_SeedImage);
	dilateFilter->SetKernel(structureElement);
	dilateFilter->SetDilateValue(1.0);
	dilateFilter->Update();

	ScalarImageTypePointer bcImage = dilateFilter->GetOutput();
	bcImage->DisconnectPipeline();

	ScalarIteratorType bciter(bcImage, bcImage->GetRequestedRegion());
	ScalarIteratorType iditer(m_voxelID, m_voxelID->GetRequestedRegion());
	while(!iditer.IsAtEnd())
	{
		if (iditer.Get()==0)
			bciter.Set(0);
        
        ++iditer;
        ++bciter;
	}

	typename BinaryContourImageFilterType::Pointer contourFilter = BinaryContourImageFilterType::New();
	contourFilter->SetInput( bcImage );
	contourFilter->SetForegroundValue( 1.0 );
	contourFilter->SetBackgroundValue( 0.0 );
	contourFilter->Update();

	ScalarIteratorType citer( contourFilter->GetOutput(), contourFilter->GetOutput()->GetRequestedRegion() );
	citer.GoToBegin();
	sditer.GoToBegin();
	while( !sditer.IsAtEnd() )
	{
		if ( citer.Get()!=0 )
			sditer.Set( 1.0 );
		else
			sditer.Set( 0.0 );

		++citer;
		++sditer;
	}

	// get foreground seeds, labeling them as 2
	std::vector<float> fgThrVec(numCc, 0.0);
	float fgThrRatio = 0.4; // top 40%
	for (int nci=0; nci<numCc; nci++)
	{	
		sort( ccSuvSet[nci].begin(), ccSuvSet[nci].end() ); // ascending order 
		int ncisize = ccSuvSet[nci].size();
		int ncind = floor(ncisize*(1.0-fgThrRatio));
		fgThrVec[nci] = ccSuvSet[nci][ncind];
	}

	nciter.GoToBegin();
	iiter.GoToBegin();
	sditer.GoToBegin();
	while( !iiter.IsAtEnd() )
	{
		int ncilab = nciter.Get();
		if (ncilab>0 && sditer.Get()==0)
		{
			float iival = iiter.Get();
			if (iival>=fgThrVec[ncilab-1] && ccSuvSet[nciter.Get()-1].size()>smallNum)
				sditer.Set(2); 
		}

		++nciter;
		++iiter;
		++sditer;
	}

}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::PerformSegmentation()
{
	if (m_logger_or_not)
		m_logger->Write("Perform segmentation ...");
	else
		std::cout << "Perform segmentation..." << std::endl;

	// create graph
	if (m_logger_or_not)
		m_logger->Write("  create graph ...");
	else
		std::cout << "  create graph ..." << std::endl;

	int wt_mat_row = m_validVxNum;
	int wt_mat_col = 1;
	for (unsigned int i=0; i<TDimension; i++)
	{
		wt_mat_col = wt_mat_col * (2*m_nei_r+1);
	}

	vnl_matrix<float> wt_mat(wt_mat_row, wt_mat_col);
	vnl_sparse_matrix<float> wt_mat_sparse(wt_mat_row, wt_mat_row);
	
	vnl_vector<float> spaDis(wt_mat_row*wt_mat_col, 0.0);
	vnl_sparse_matrix<float> spaDis_mat_sp(wt_mat_row, wt_mat_row);

	vnl_vector<float> feaDis(wt_mat_row*wt_mat_col, 0.0);
	std::vector<vnl_vector<float> > feaDisVec(TInputImageNum, feaDis);
	vnl_sparse_matrix<float> dis_mat_sparse(wt_mat_row, wt_mat_row);
	std::vector<vnl_sparse_matrix<float> > feaDis_mat_sp_vec(TInputImageNum, dis_mat_sparse);

	VectorNeighborhoodRadiusType vec_nei_r;
	ScalarNeighborhoodRadiusType sca_nei_r;
	for (unsigned int i=0; i<TDimension; i++)
	{
		vec_nei_r[i] = m_nei_r;
		sca_nei_r[i] = m_nei_r;
	}

	VectorNeighborhoodIteratorType vecImgIter(vec_nei_r, m_InputImageVector, m_InputImageVector->GetRequestedRegion());
	ScalarNeighborhoodIteratorType voxelIdIter(sca_nei_r, m_voxelID, m_voxelID->GetRequestedRegion());

	typename VectorImageType::PixelType ftemp;
	typename VectorImageType::IndexType stemp;
	int cenId = -1;
	int neiId = -1;
	float numEdge = 0;

	// // get distance in feature/spatial space
	vecImgIter.GoToBegin();
	voxelIdIter.GoToBegin();
	while( !vecImgIter.IsAtEnd() )
	{
		cenId = voxelIdIter.GetCenterPixel() - 1;
		if (cenId>=0)
		{
			for (int i=0; i<wt_mat_col; i++)
			{
				neiId = (int)(voxelIdIter.GetPixel(i)-1);
				if (neiId>=0)
				{
					numEdge++;

					ftemp = vecImgIter.GetCenterPixel() - vecImgIter.GetPixel(i);
					for (int mi=0; mi<TInputImageNum; mi++)
					{
						feaDisVec[mi][cenId*wt_mat_col+i] = ftemp[mi] * ftemp[mi];
						feaDis_mat_sp_vec[mi](cenId,neiId) = ftemp[mi] * ftemp[mi];
					}

					float s2 = 0;
					for (int ind=0; ind<TDimension; ind++)
					{
						stemp[ind] = vecImgIter.GetIndex()[ind] - vecImgIter.GetIndex(i)[ind];
						s2 = s2 + (stemp[ind]*m_ImageVoxelSize[ind])*(stemp[ind]*m_ImageVoxelSize[ind]);
					}						
					spaDis[cenId*wt_mat_col+i] = s2;
					spaDis_mat_sp(cenId,neiId) = s2;
				}
			}
		}
		++vecImgIter;
		++voxelIdIter;
	}
	// // get similarity sigma
	//cout << m_validVxNum << "	" << numEdge << endl;
	//cout << wt_mat_row << "    " << wt_mat_col << endl;

	std::vector<float> feaSig(TInputImageNum, 0.0);
	for (unsigned int ti=0; ti<TInputImageNum; ti++)
		feaSig[ti] = feaDisVec[ti].sum() / numEdge;
	//cout << feaSig[0] << "	" << feaSig[1] << endl;

	float spaSig = spaDis.sum() / numEdge;
	//cout << spaSig << endl;

	voxelIdIter.GoToBegin();
	cenId = -1;
	neiId = -1;
	while( !voxelIdIter.IsAtEnd() )
	{
		cenId = voxelIdIter.GetCenterPixel() - 1;
		if ( cenId>=0 )
		{
			for (int i=0; i<wt_mat_col; i++)
			{
				neiId = (int)(voxelIdIter.GetPixel(i)-1);
				if (neiId>=0)
				{
					float spaWt = exp(-spaDis_mat_sp(cenId,neiId)/spaSig);
					float feaWt = 0.0;
					for (unsigned ti=0; ti<TInputImageNum; ti++)
						feaWt += exp(-feaDis_mat_sp_vec[ti](cenId,neiId)/feaSig[ti]) * m_modailityWt[ti];

					wt_mat_sparse(cenId,neiId) = feaWt * spaWt;
				}
			}
			wt_mat_sparse(cenId,cenId) = 0.0;
		}
		++voxelIdIter;
	}
	// // similarity normalization
	voxelIdIter.GoToBegin();
	cenId = -1;
	neiId = -1;
	while( !voxelIdIter.IsAtEnd() )
	{
		cenId = voxelIdIter.GetCenterPixel() - 1;
		if ( cenId>=0 )
		{
			for (int i=0; i<wt_mat_col; i++)
			{
				neiId = (int)(voxelIdIter.GetPixel(i)-1);
				if (neiId>=0)
					wt_mat(cenId,i) = wt_mat_sparse(cenId,neiId);
				else
					wt_mat(cenId,i) = 0.0;
			}
		}
		++voxelIdIter;
	}

	vnl_vector<float> invSqrtRowSum(wt_mat.rows(),0.0);
	for (unsigned int i=0; i<wt_mat.rows(); i++)
	{
		float wt_i_sum = wt_mat.get_row(i).sum();
		if (wt_i_sum>0)
		{
			invSqrtRowSum[i] = 1.0 / vcl_sqrt(wt_i_sum);
		}
		//cout << i << ". " << invSqrtRowSum[i] << endl;
	}
 
 	voxelIdIter.GoToBegin();
	cenId = -1;
	neiId = -1;
	while( !voxelIdIter.IsAtEnd() )
	{
		cenId = voxelIdIter.GetCenterPixel() - 1;
		if ( cenId>=0 )
		{
			//cout << cenId << ": ";
			for (int i=0; i<wt_mat_col; i++)
			{
				neiId = (int)(voxelIdIter.GetPixel(i)-1);
				if (neiId>=0)
					wt_mat(cenId,i) = wt_mat(cenId,i) * invSqrtRowSum[cenId] * invSqrtRowSum[neiId];
				//cout << wt_mat(cenId,i) << ", ";
			}
			//cout << endl;
		}
		++voxelIdIter;
	}

	// label initialization
	if (m_logger_or_not)
		m_logger->Write("  label initialization ...");
	else
		std::cout << "  label initialization ..." << std::endl;
	
	LabValVectorType tmpLabVec(0.0);

	LabValImageTypePointer iniLabImg = LabValImageType::New();
	iniLabImg->SetOrigin( m_MfImage[0]->GetOrigin() );
	iniLabImg->SetSpacing( m_MfImage[0]->GetSpacing() );
	iniLabImg->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	iniLabImg->SetDirection( m_MfImage[0]->GetDirection() );
	iniLabImg->Allocate(); 
	iniLabImg->FillBuffer( tmpLabVec );

	LabValIteratorType labValIter(iniLabImg, iniLabImg->GetRequestedRegion());
	ScalarIteratorType iniLabIter(m_SeedImage, m_SeedImage->GetRequestedRegion());
	ScalarIteratorType vxlIdIter(m_voxelID, m_voxelID->GetRequestedRegion());

	labValIter.GoToBegin();
	iniLabIter.GoToBegin();
	vxlIdIter.GoToBegin();
	while( !vxlIdIter.IsAtEnd() )
	{
		if (vxlIdIter.Get()>0)
		{
			int lab = (int)iniLabIter.Get();
			
			LabValVectorType tmpLab(0.0);
			if (lab!=0)
			{
				tmpLab[lab-1] = 1.0;
			}
			labValIter.Set(tmpLab);
		}

		++labValIter;
		++iniLabIter;
		++vxlIdIter;
	}

	// label propagation
	if (m_logger_or_not)
		m_logger->Write("  label propagation ...");
	else
		std::cout << "  label propagation ..." << std::endl;
	
	LabValImageTypePointer tempLabValImg = LabValImageType::New();
	tempLabValImg->SetOrigin( m_MfImage[0]->GetOrigin() );
	tempLabValImg->SetSpacing( m_MfImage[0]->GetSpacing() );
	tempLabValImg->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	tempLabValImg->SetDirection( m_MfImage[0]->GetDirection() );
	tempLabValImg->Allocate(); 
	tempLabValImg->FillBuffer( tmpLabVec );

	LabValImageTypePointer labValImg = LabValImageType::New();
	labValImg->SetOrigin( m_MfImage[0]->GetOrigin() );
	labValImg->SetSpacing( m_MfImage[0]->GetSpacing() );
	labValImg->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	labValImg->SetDirection( m_MfImage[0]->GetDirection() );
	labValImg->Allocate(); 
	labValImg->FillBuffer( tmpLabVec );

	LabValNeighborhoodRadiusType labVal_r;
	for (unsigned int i=0; i<TDimension; i++)
	{
		labVal_r[i] = m_nei_r;
	}

	LabValNeighborhoodIteratorType initer(labVal_r, iniLabImg, iniLabImg->GetRequestedRegion());
	LabValNeighborhoodIteratorType tmpiter(labVal_r, tempLabValImg, tempLabValImg->GetRequestedRegion());
	LabValNeighborhoodIteratorType valiter(labVal_r, labValImg, labValImg->GetRequestedRegion());	
	//ScalarNeighborhoodIteratorType voxelIdIter(labVal_r, m_voxelID, m_voxelID->GetRequestedRegion());

	tmpiter.GoToBegin();
	initer.GoToBegin();
	while( !initer.IsAtEnd() )
	{
		tmpiter.SetCenterPixel(initer.GetCenterPixel());

		++tmpiter;
		++initer;
	}

	int iter = 0;
	while( iter<m_iterNum )
	{
		initer.GoToBegin();
		tmpiter.GoToBegin();
		valiter.GoToBegin();
		voxelIdIter.GoToBegin();

		while( !voxelIdIter.IsAtEnd() )
		{
			int cid = voxelIdIter.GetCenterPixel() - 1;
			if (cid>=0)
			{
				LabValVectorType prop_tmp(0.0);
				for (int c=0; c<wt_mat_col; c++)
				{
					prop_tmp = prop_tmp + tmpiter.GetPixel(c) * wt_mat(cid,c) * m_alpha;
				}
				prop_tmp = prop_tmp + initer.GetCenterPixel() * (1-m_alpha);

				valiter.SetCenterPixel(prop_tmp);
			}

			++initer;
			++tmpiter;
			++valiter;
			++voxelIdIter;
		}

		tmpiter.GoToBegin();
		valiter.GoToBegin();
		while(!tmpiter.IsAtEnd())
		{
			tmpiter.SetCenterPixel(valiter.GetCenterPixel());
			++tmpiter;
			++valiter;
		}

		iter++;
		if (iter % 50 == 0)
		{
			if (m_logger_or_not)
				m_logger->Write("    iteration " + std::to_string(iter) + " in " + std::to_string(m_iterNum) + " completed");
			else
				std::cout << "    iteration " << iter << " in " << m_iterNum << " completed" << std::endl;
		}
	}

	// get segmentation label
	if (m_logger_or_not)
		m_logger->Write("  get segmentation label ...");
	else
		std::cout << "  get segmentation label ..." << std::endl;

	m_segImage = ScalarLabelImageType::New();
	m_segImage->SetOrigin( m_MfImage[0]->GetOrigin() );
	m_segImage->SetSpacing( m_MfImage[0]->GetSpacing() );
	m_segImage->SetRegions( m_MfImage[0]->GetRequestedRegion() );
	m_segImage->SetDirection( m_MfImage[0]->GetDirection() );
	m_segImage->Allocate(); 
	m_segImage->FillBuffer( 0 );

	ScalarLabelImageIteratorType segiter(m_segImage, m_segImage->GetRequestedRegion());
	ScalarIteratorType svaliter(m_OutputImage, m_OutputImage->GetRequestedRegion());

	svaliter.GoToBegin();
	segiter.GoToBegin();
	valiter.GoToBegin();
	voxelIdIter.GoToBegin();
	while( !voxelIdIter.IsAtEnd() )
	{
		if (voxelIdIter.GetCenterPixel()>0)
		{
			TPixelType maxVal = 0;
			unsigned int maxInd = 0;
			for (unsigned int ci=0; ci<2; ci++)
			{
				if (valiter.GetCenterPixel()[ci]>maxVal)
				{
					maxVal = valiter.GetCenterPixel()[ci];
					maxInd = ci;
				}
			}
			if (maxVal>0)
				segiter.Set(maxInd+1);
            if (maxVal==0)
                segiter.Set(1);

			svaliter.Set(valiter.GetCenterPixel()[1]);
		}

		++svaliter;
		++segiter;
		++valiter;
		++voxelIdIter;
	}

    // relabel image
    segiter.GoToBegin();
    voxelIdIter.GoToBegin();
    while( !voxelIdIter.IsAtEnd() )
    {
        if (segiter.Get()==2)
            segiter.Set(1);
        else if (segiter.Get()==1)
            segiter.Set(2);
        
        if (voxelIdIter.GetCenterPixel()==-1)
        {
            segiter.Set(3);
        }
        
        ++segiter;
        ++voxelIdIter;
    }
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_Nodule< TPixelType, TDimension, TInputImageNum >::WriteLabel(int seedAv)
{
	if (m_logger_or_not)
		m_logger->Write("Write label");
	else
		std::cout << "Write label" << std::endl;

    if (seedAv==0)
    {
        std::string outName = m_outputBaseName + "_seeds.nii.gz";

        typename ScalarImageWriterType::Pointer labelWriter = ScalarImageWriterType::New();
        labelWriter->SetFileName( outName );
        labelWriter->SetInput( m_SeedImage );
        labelWriter->Update();
    }
    
	//string outValName = m_outputBaseName + "_lab_val.nii.gz";

	//labelWriter->SetFileName( outValName );
	//labelWriter->SetInput( m_OutputImage );
	//labelWriter->Update();

	std::string labelingName = m_outputBaseName + "_segmentation.nii.gz";

	typename ScalarLabelImageWriterType::Pointer labelingWriter = ScalarLabelImageWriterType::New();
	labelingWriter->SetFileName( labelingName );
	labelingWriter->SetInput( m_segImage );
	labelingWriter->Update();
}


