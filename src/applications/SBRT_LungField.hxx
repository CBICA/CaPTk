#pragma once

#include <sstream>
#include "SBRT_LungField.h"

using namespace std;
using namespace cv;

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
SBRT_LungField< TPixelType, TDimension, TInputImageNum >::SBRT_LungField()
{
  m_maskAv = 0;
	m_compactness = 1;
	m_iterNum = 100;

	m_logger = NULL;
	m_logger_or_not = false;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
SBRT_LungField< TPixelType, TDimension, TInputImageNum >::~SBRT_LungField()
{

}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::ReadMask( string maskName )
{
	typename ScalarImageReaderType::Pointer MaskReader = ScalarImageReaderType::New();
	MaskReader->SetFileName(maskName);
	MaskReader->Update();

  ScalarIteratorType Riter(MaskReader->GetOutput(), MaskReader->GetOutput()->GetRequestedRegion());
  float tmpSum = 0.0;
  Riter.GoToBegin();
  while (!Riter.IsAtEnd())
  {
    if (Riter.Get() > 0)
    {
      tmpSum += 1;
    }

    ++Riter;
  }

  if (tmpSum == 0)
    m_maskAv = 0;
  else
  {
    m_Mask = ScalarImageType::New();
    m_Mask->SetOrigin(MaskReader->GetOutput()->GetOrigin());
    m_Mask->SetSpacing(MaskReader->GetOutput()->GetSpacing());
    m_Mask->SetRegions(MaskReader->GetOutput()->GetRequestedRegion());
    m_Mask->SetDirection(MaskReader->GetOutput()->GetDirection());
    m_Mask->Allocate();
    m_Mask->FillBuffer(0.0);

    ScalarIteratorType Miter(m_Mask, m_Mask->GetRequestedRegion());
    Riter.GoToBegin();
    Miter.GoToBegin();
    while (!Riter.IsAtEnd())
    {
      if (Riter.Get() > 0)
      {
        Miter.Set(1.0);
      }

      ++Riter;
      ++Miter;
    }
  }
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
int SBRT_LungField< TPixelType, TDimension, TInputImageNum >::Initialization( vector<string> &InputImageName )
{
	if (m_logger_or_not)
		m_logger->Write("Read images ...");
	else
		std::cout << "Read images... " << std::endl;	
	
	typename ScalarImageReaderType::Pointer InputReader[TInputImageNum];
	for ( unsigned int i=0;i<TInputImageNum;i++ )
	{
		InputReader[i] = ScalarImageReaderType::New();
		InputReader[i]->SetFileName( InputImageName[i] );

		try
		{
			InputReader[i]->Update();
		}
		catch ( itk::ExceptionObject &err )
		{
			std::cerr << " ExceptionObject caught !" << std::endl;
			std::cerr << " " << err << std::endl;
			return EXIT_FAILURE;
		}
	}

	m_InputImageVector = VectorImageType::New();
	m_InputImageVector->SetOrigin( InputReader[0]->GetOutput()->GetOrigin() );
	m_InputImageVector->SetSpacing( InputReader[0]->GetOutput()->GetSpacing() );
	m_InputImageVector->SetRegions( InputReader[0]->GetOutput()->GetRequestedRegion() );
	m_InputImageVector->SetDirection( InputReader[0]->GetOutput()->GetDirection() );
	m_InputImageVector->Allocate();
	VectorType temp(-1000); 
	m_InputImageVector->FillBuffer( temp );

	VectorIteratorType Viter( m_InputImageVector,m_InputImageVector->GetRequestedRegion() );
	for ( unsigned int i=0; i<TInputImageNum; i++ )
	{
		ScalarIteratorType Siter( InputReader[i]->GetOutput(),InputReader[i]->GetOutput()->GetRequestedRegion() );
		Siter.GoToBegin();
		Viter.GoToBegin();
		while( !Viter.IsAtEnd() )
		{
			temp = Viter.Get();
			temp[i] = Siter.Get();
			Viter.Set(temp);				
			
			++Viter;
			++Siter;
		}
	}

	// get foreground mask
	if (m_maskAv ==0)
	{
		m_Mask = ScalarImageType::New();
		m_Mask->SetOrigin( InputReader[0]->GetOutput()->GetOrigin() );
		m_Mask->SetSpacing( InputReader[0]->GetOutput()->GetSpacing() );
		m_Mask->SetRegions( InputReader[0]->GetOutput()->GetRequestedRegion() );
		m_Mask->SetDirection( InputReader[0]->GetOutput()->GetDirection() );
		m_Mask->Allocate(); 
		m_Mask->FillBuffer( 0 );

		GetImageMask();		
	}

	// mask input images
	VectorType zeroVec(0.0);
	ScalarIteratorType Miter( m_Mask, m_Mask->GetRequestedRegion() );
	
	Viter.GoToBegin();
	Miter.GoToBegin();
	while( !Viter.IsAtEnd() )
	{
		if ( Miter.Get()==0 )
		{
			Viter.Set(zeroVec);				
		}
		
		++Viter;
		++Miter;
	}

	// calculate image gradient
	m_GraMagnitude = ScalarImageType::New();
	m_GraMagnitude->SetOrigin( InputReader[0]->GetOutput()->GetOrigin() );
	m_GraMagnitude->SetSpacing( InputReader[0]->GetOutput()->GetSpacing() );
	m_GraMagnitude->SetRegions( InputReader[0]->GetOutput()->GetRequestedRegion() );
	m_GraMagnitude->SetDirection( InputReader[0]->GetOutput()->GetDirection() );
	m_GraMagnitude->Allocate(); 
	m_GraMagnitude->FillBuffer( -1 );

	if ( 1==TInputImageNum )
	{
		typename ScalarGraMagFilterType::Pointer gmFilter = ScalarGraMagFilterType::New();
		gmFilter->SetInput(InputReader[0]->GetOutput());

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
	m_OutputImage->SetOrigin( InputReader[0]->GetOutput()->GetOrigin() );
	m_OutputImage->SetSpacing( InputReader[0]->GetOutput()->GetSpacing() );
	m_OutputImage->SetRegions( InputReader[0]->GetOutput()->GetRequestedRegion() );
	m_OutputImage->SetDirection( InputReader[0]->GetOutput()->GetDirection() );
	m_OutputImage->Allocate(); 
	m_OutputImage->FillBuffer( -1 );
	// initialize valid voxel ID image
	m_voxelID = ScalarImageType::New();
	m_voxelID->SetOrigin( InputReader[0]->GetOutput()->GetOrigin() );
	m_voxelID->SetSpacing( InputReader[0]->GetOutput()->GetSpacing() );
	m_voxelID->SetRegions( InputReader[0]->GetOutput()->GetRequestedRegion() );
	m_voxelID->SetDirection( InputReader[0]->GetOutput()->GetDirection() );
	m_voxelID->Allocate(); 
	m_voxelID->FillBuffer( 0 );

	//ScalarIteratorType Siter( m_OutputImage, m_OutputImage->GetRequestedRegion() );
	ScalarIteratorType Iiter( m_voxelID, m_voxelID->GetRequestedRegion() );
	//Siter.GoToBegin();
	Iiter.GoToBegin();
	Viter.GoToBegin();
	m_validVxNum = 0;
	while ( !Viter.IsAtEnd() )
	{		
		if ( Viter.Get()[0]>0 )
		{
			m_validVxNum++;
			Iiter.Set( m_validVxNum );
		}

		++Iiter;
		//++Siter;
		++Viter;
	}
	if (m_logger_or_not)
		m_logger->Write("	valid voxel number: " + to_string(m_validVxNum));
	else
		std::cout << "	valid voxel number: " << m_validVxNum << std::endl;

	for( unsigned int j=0;j<TDimension;j++ )
	{
		m_ImageVoxelSize.push_back( m_InputImageVector->GetSpacing()[j] );
	}

	m_ImageSize = m_InputImageVector->GetRequestedRegion().GetSize();

	//cout << "Read images successfully " << endl;
    
    return 0;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::GetImageMask()
{
	float petThr = 0.05;
	float ctThr = -950;

	VectorIteratorType Viter( m_InputImageVector,m_InputImageVector->GetRequestedRegion() );
	ScalarIteratorType Miter( m_Mask, m_Mask->GetRequestedRegion() );

	Viter.GoToBegin();
	Miter.GoToBegin();
	while(!Viter.IsAtEnd())
	{
		if ( Viter.Get()[0]>petThr && Viter.Get()[1]>ctThr )
			Miter.Set(1.0);
		
		++Viter;
		++Miter;
	}

	typename BinaryFunctionType::Pointer binaryFunction = BinaryFunctionType::New();
	binaryFunction->SetInputImage(m_Mask);
	binaryFunction->ThresholdAbove(1.0);

	float regionLabel = 2.0;
	float maxNum = 0;
	float maxNumLabel = 0;
	Miter.GoToBegin();
	while (!Miter.IsAtEnd())
	{
		if (Miter.Get() == 1.0)
		{
			FloodFilledIterator floodIterator(m_Mask,binaryFunction,Miter.GetIndex());
			floodIterator.GoToBegin();
			float counter = 0;
			while(!floodIterator.IsAtEnd())
			{
				floodIterator.Set(regionLabel);
				++floodIterator;
				counter++;
			}
	
			if (maxNum < counter)
			{
				maxNum = counter;
				maxNumLabel = regionLabel;
			}
			regionLabel++;
		}
		++Miter;
	}
	
	Miter.GoToBegin();
	while (!Miter.IsAtEnd())
	{
		if (Miter.Get() != maxNumLabel)
			Miter.Set(0.0);
		else
			Miter.Set(1.0);

		++Miter;
	}
	
	
	if (m_logger_or_not)
		m_logger->Write("Erode and dilate ...");
	else
		std::cout << "Erode and dilate ..." << std::endl;
	
	typedef itk::BinaryBallStructuringElement<TPixelType,TDimension> StructureElementType;
	StructureElementType structureElement;
	structureElement.SetRadius(3);
	structureElement.CreateStructuringElement();
	
	typedef itk::BinaryErodeImageFilter<ScalarImageType,ScalarImageType,StructureElementType> BinaryErodeImageFilterType;
	typename BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
	erodeFilter->SetInput(m_Mask);
	erodeFilter->SetKernel(structureElement);
	erodeFilter->SetErodeValue(1.0);
	erodeFilter->Update();

	structureElement.SetRadius(1);
	structureElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter<ScalarImageType,ScalarImageType,StructureElementType> BinaryDilateImageFilterType;
	typename BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput(erodeFilter->GetOutput());
	dilateFilter->SetKernel(structureElement);
	dilateFilter->SetDilateValue(1.0);
	dilateFilter->Update();

	ScalarIteratorType DilateResultIterator(dilateFilter->GetOutput(), dilateFilter->GetOutput()->GetRequestedRegion());
	DilateResultIterator.GoToBegin();
	Miter.GoToBegin();
	while(!Miter.IsAtEnd())
	{
		if (DilateResultIterator.Get()==1.0)
			Miter.Set(1.0);
		else
			Miter.Set(0.0);

		++Miter;
		++DilateResultIterator;
	}
	
	// pick the largerst one
	maxNum = 0;
	maxNumLabel = 0;
	Miter.GoToBegin();
	while (!Miter.IsAtEnd())
	{
		if (Miter.Get() == 1.0)
		{
			FloodFilledIterator floodIterator(m_Mask, binaryFunction, Miter.GetIndex());
			floodIterator.GoToBegin();
			float counter = 0;
			while (!floodIterator.IsAtEnd())
			{
				floodIterator.Set(regionLabel);
				++floodIterator;
				counter++;
			}

			if (maxNum < counter)
			{
				maxNum = counter;
				maxNumLabel = regionLabel;
			}
			regionLabel++;
		}
		++Miter;
	}

	Miter.GoToBegin();
	while (!Miter.IsAtEnd())
	{
		if (Miter.Get() != maxNumLabel)
			Miter.Set(0.0);
		else
			Miter.Set(1.0);

		++Miter;
	}
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::SetParameters(int maskAv, float svxSize, float compactness, float minSize, int iter, int ptSz )
{
  m_maskAv = maskAv;
	m_svxSize = svxSize;
	m_compactness = compactness;
	m_minSize = minSize;
	m_iterNum = iter;
	m_ptSz = ptSz;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::SetLogger(string logName)
{
	m_logger = new cbica::Logging;
	m_logger->UseNewFile(logName);
	m_logger_or_not = true;
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::GenerateSeeds()
{
	if (m_logger_or_not)
		m_logger->Write("Generate seeds ...");
	else
		std::cout << "Generate seeds ..." << std::endl;
	
	m_step = pow( m_svxSize,1.0f/TDimension ) + 0.5;

	if ( TDimension==3 )
	{
		int xstrips = ( 0.5+float(m_ImageSize[0])/float(m_step) );
		int ystrips = ( 0.5+float(m_ImageSize[1])/float(m_step) );
		int zstrips = ( 0.5+float(m_ImageSize[2])/float(m_step) );

		int xerr = m_ImageSize[0] - m_step*xstrips;
		if(xerr < 0){ xstrips--; xerr = m_ImageSize[0] - m_step*xstrips;}
		int yerr = m_ImageSize[1] - m_step*ystrips;
		if(yerr < 0){ ystrips--; yerr = m_ImageSize[1] - m_step*ystrips;}
		int zerr = m_ImageSize[2] - m_step*zstrips;
		if(zerr < 0){ zstrips--; zerr = m_ImageSize[2] - m_step*zstrips;}

		double xerrperstrip = double(xerr)/double(xstrips);
		double yerrperstrip = double(yerr)/double(ystrips);
		double zerrperstrip = double(zerr)/double(zstrips);

		int xoff = m_step/2;
		int yoff = m_step/2;
		int zoff = m_step/2;

		
		for ( int z=0;z<zstrips;z++ )
		{
			int ze = z*zerrperstrip;
			int zp = z*m_step+zoff+ze;
			for ( int y=0;y<ystrips;y++ )
			{
				int ye = y*yerrperstrip;
				int yp = y*m_step+yoff+ye;
				for ( int x=0;x<xstrips;x++ )
				{
					int xe = x*xerrperstrip;
					int xp = x*m_step+xoff+xe;

					VectorIndexType curInd;
					curInd[0] = xp;
					curInd[1] = yp;
					curInd[2] = zp;
					//cout << "  " << xp << "  " << yp << "  " << zp << endl;
					if ( m_voxelID->GetPixel(curInd)>0 )
					{
						// adjust the seed according to image gradient
						float curGm = m_GraMagnitude->GetPixel(curInd);
						for (int nz=zp-1;nz<zp+1;nz++)
							for (int ny=yp-1;ny<yp+1;ny++)
								for (int nx=xp-1;nx<xp+1;nx++)
								{
									if ( (nx>=0 && nx<m_ImageSize[0]) && (ny>=0 && ny<m_ImageSize[1]) && (nz>=0 && nz<m_ImageSize[2]) )
									{
										VectorIndexType curNind;
										curNind[0] = nx;
										curNind[1] = ny;
										curNind[2] = nz;
										if ( m_voxelID->GetPixel(curNind)>0 )
										{
											float curNgm = m_GraMagnitude->GetPixel(curNind);
											if (curNgm<curGm)
											{
												curInd[0] = nx;
												curInd[1] = ny;
												curInd[2] = nz;
												curGm = curNgm;
											}
										}
									}
								}

						// assign the adjusted seed
						int cInd = m_voxelID->GetPixel(curInd) - 1;
						vnl_vector<float> curImgVal = m_ImageFea.get_row(cInd);
						for ( int ii=0;ii<m_feaNum;ii++ )
						{
							m_SeedCentroids[ii].push_back(curImgVal[ii]);
						}

						m_SeedCentroids[m_feaNum].push_back(curInd[0]);
						m_SeedCentroids[m_feaNum+1].push_back(curInd[1]);
						m_SeedCentroids[m_feaNum+2].push_back(curInd[2]);
					}
				}
			}
		}
	}
	else
	{
		if (m_logger_or_not)
			m_logger->Write("  support 3D image only now !");
		else
			std::cout << "  support 3D image only now !" << std::endl;
	}
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::FeaExtraction()
{
	if (m_logger_or_not)
		m_logger->Write("Feature extraction ...");
	else
		std::cout << "Feature extraction..." << std::endl;
	
	int NeighNum = (int)pow((double)(2*m_ptSz+1),(double)TDimension);
	m_feaNum = NeighNum * TInputImageNum;
	m_ImageFea.set_size(m_validVxNum,m_feaNum);
	m_ImageFea.fill(0.0f);

	VectorNeighborhoodRadiusType vectorRadius;
	for (int di=0;di<TDimension;di++)
		vectorRadius[di] = m_ptSz;
	
	VectorNeighborhoodIteratorType VecNeighIter( vectorRadius,m_InputImageVector,m_InputImageVector->GetRequestedRegion() );
	ScalarIteratorType Iiter( m_voxelID,m_voxelID->GetRequestedRegion() );
	VecNeighIter.GoToBegin();
	Iiter.GoToBegin();
	while( !Iiter.IsAtEnd() )
	{
		if (Iiter.Get()>0)
		{
			int vInd = Iiter.Get()-1; // index from 0
			int fInd = 0;
			for ( int i=0;i<NeighNum;i++ )
			{
				for( int k=0;k<TInputImageNum;k++ )
				{
					m_ImageFea(vInd,fInd) = VecNeighIter.GetPixel(i)[k];
					fInd++;
				}
			}
		}
		++Iiter;
		++VecNeighIter;
	}

	m_ImageFea.normalize_columns();
    m_ImageFea *= 10.0; 
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::PerformSLIC( float compactness )
{
	if (m_logger_or_not)
		m_logger->Write("Perform SLIC ...");
	else
		std::cout << "Perform SLIC..." << std::endl;

	int numK = m_SeedCentroids[0].size();
	
	// for 3D only now
	//------------------
	int offset = m_step;
	//------------------
	vector<float> clusterSize(numK,0);
	vector<float> inv(numK,0);
	vector<float> vecFlt(numK,0);
	vector<float> distvec(m_validVxNum,FLT_MAX);
	vector<float> maxFeaDis(numK,10);
	float invwt = compactness/pow(m_step,2.0);

	int x1,x2,y1,y2,z1,z2;
	vnl_vector<float> tmpVal;
	float distfea;
	float distxyz;
	float energy;
	float oldEnergy;
	float startEnergy;
	for ( int iter=0;iter<m_iterNum;iter++ ) // how many iterations
	{
		//cout << "  iteration: " << iter+1 << endl;
	//	std::stringstream ss;
	//	ss << "  iteration: " << iter + 1 << endl;
	//		if (m_logger_or_not)
	//	m_logger->Write(ss.str());
	//else
	//	cout << ss.str() << endl;
		vector<float> tmpMaxFeaDis(numK,0);
		energy = 0;

		distvec.assign( m_validVxNum,FLT_MAX );
		//cout << "    update membership" << endl;
		for ( int n=0;n<numK;n++ )
		{
			y1 = max( 0,int(m_SeedCentroids[m_feaNum+1][n]-offset) );
			y2 = min( int(m_ImageSize[1]),int(m_SeedCentroids[m_feaNum+1][n]+offset) );
			x1 = max( 0,int(m_SeedCentroids[m_feaNum][n]-offset) );
			x2 = min( int(m_ImageSize[0]),int(m_SeedCentroids[m_feaNum][n]+offset) );
			z1 = max( 0,int(m_SeedCentroids[m_feaNum+2][n]-offset) );
			z2 = min( int(m_ImageSize[2]),int(m_SeedCentroids[m_feaNum+2][n]+offset) );

			for ( int z=z1;z<z2;z++ )
			{
				for ( int y=y1;y<y2;y++ )
				{
					for ( int x=x1;x<x2;x++ )
					{
						VectorIndexType curInd;
						curInd[0] = x;
						curInd[1] = y;
						curInd[2] = z;

						int vxId = m_voxelID->GetPixel(curInd);
						if ( vxId>0 )
						{
							int vInd = vxId - 1;
							tmpVal = m_ImageFea.get_row(vInd);
							
							distfea = 0;
							for ( int fi=0;fi<m_feaNum;fi++ )
							{
								distfea += (tmpVal(fi)-m_SeedCentroids[fi][n])*(tmpVal(fi)-m_SeedCentroids[fi][n]);
							}

							distxyz = (x-m_SeedCentroids[m_feaNum][n])*(x-m_SeedCentroids[m_feaNum][n])+
								      (y-m_SeedCentroids[1+m_feaNum][n])*(y-m_SeedCentroids[1+m_feaNum][n])+
									  (z-m_SeedCentroids[2+m_feaNum][n])*(z-m_SeedCentroids[2+m_feaNum][n]);
							
							float dist = distfea/maxFeaDis[n] + distxyz*invwt;
							
							if ( dist<distvec[vxId-1] )
							{
								distvec[vxId-1] = dist;
								m_OutputImage->SetPixel( curInd,n );

								if ( tmpMaxFeaDis[n]<distfea )
									tmpMaxFeaDis[n] = distfea;
							}
						}
					}
				}
			}
		}
		for ( int n=0;n<numK;n++ )
			maxFeaDis[n] = tmpMaxFeaDis[n];

		for ( unsigned int di=0;di<distvec.size();di++ )
			if ( distvec[di]!=FLT_MAX )
				energy += distvec[di];
		if ( iter==0 )
			startEnergy = energy;
		else
		{
			//if ( abs(oldEnergy-energy)<1e-6*abs(startEnergy-energy) )
			if ( abs(oldEnergy-energy)<1e-6*oldEnergy )
				break;
		}
		oldEnergy = energy;

		// recalculate the centroid and store the seed values
		//cout << "    update seed centroids..." << endl;
		vector<vector<float> > sigma(TDimension+m_feaNum,vecFlt);
		clusterSize.assign( numK,0 );
		for ( int z=0;z<m_ImageSize[2];z++ )
		{
			for ( int y=0;y<m_ImageSize[1];y++ )
			{
				for ( int x=0;x<m_ImageSize[0];x++ )
				{
					VectorIndexType curInd;
					curInd[0] = x;
					curInd[1] = y;
					curInd[2] = z;

					int vId = m_voxelID->GetPixel(curInd);
					if ( vId>0 )
					{
						vnl_vector<float> curVal = m_ImageFea.get_row(vId-1);
						int curLab = m_OutputImage->GetPixel(curInd);
						
						if ( curLab>=0 )
						{
							for ( int fi=0;fi<m_feaNum;fi++ )
							{
								sigma[fi][curLab] += curVal(fi);
							}
							sigma[m_feaNum][curLab] += x;
							sigma[m_feaNum+1][curLab] += y;
							sigma[m_feaNum+2][curLab] += z;

							clusterSize[curLab] += 1;
						}
						
					}
				}
			}
		}
		for ( int k=0;k<numK;k++ )
		{
			if ( clusterSize[k]<=0 ) clusterSize[k] = 1;
			inv[k] = 1/clusterSize[k];
		}
		for ( int k=0;k<numK;k++ )
		{
			for ( int fi=0;fi<m_feaNum;fi++ )
			{
				m_SeedCentroids[fi][k] = sigma[fi][k] * inv[k];
			}
			m_SeedCentroids[m_feaNum][k] = sigma[m_feaNum][k] * inv[k];
			m_SeedCentroids[m_feaNum+1][k] = sigma[m_feaNum+1][k] * inv[k];
			m_SeedCentroids[m_feaNum+2][k] = sigma[m_feaNum+2][k] * inv[k];
		}
	}

}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::EnforceSupervoxelLabelConnectivity()
{
	if (m_logger_or_not)
		m_logger->Write("Enforce label connectivity ...");
	else
		std::cout << "Enforce label connectivity..." << std::endl;

	const int dx10[10] = {-1,  0,  1,  0, -1,  1,  1, -1,  0, 0};
	const int dy10[10] = { 0, -1,  0,  1, -1, -1,  1,  1,  0, 0};
	const int dz10[10] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, 1};

	int adjlabel = 0; //adjacent label
	vector<int> xvec;
	vector<int> yvec;
	vector<int> zvec;
	
	m_labelingImage = ScalarLabelImageType::New();
	m_labelingImage->SetOrigin( m_OutputImage->GetOrigin() );
	m_labelingImage->SetSpacing( m_OutputImage->GetSpacing() );
	m_labelingImage->SetRegions( m_OutputImage->GetRequestedRegion() );
	m_labelingImage->SetDirection( m_OutputImage->GetDirection() );
	m_labelingImage->Allocate(); 
	m_labelingImage->FillBuffer( -1 );
	//------------------
	// labeling
	//------------------
	int lab = 0;
	for( int d=0; d<m_ImageSize[2]; d++ )
	{
		for( int h=0; h<m_ImageSize[1]; h++ )
		{
			for( int w=0; w<m_ImageSize[0]; w++ )
			{
				VectorIndexType curInd;
				curInd[0] = w;
				curInd[1] = h;
				curInd[2] = d;

				if ( m_voxelID->GetPixel(curInd)>0 )
				{
					if( m_labelingImage->GetPixel(curInd)<0 )
					{
						m_labelingImage->SetPixel( curInd,lab );
						//-------------------------------------------------------
						// Quickly find an adjacent label for use later if needed
						//-------------------------------------------------------
						for( int n=0; n<10; n++ )
						{
							int x = w + dx10[n];
							int y = h + dy10[n];
							int z = d + dz10[n];
							if( (x>=0 && x<m_ImageSize[0]) && (y>=0 && y<m_ImageSize[1]) && (z>=0 && z<m_ImageSize[2]) )
							{
								VectorIndexType curNind;
								curNind[0] = x;
								curNind[1] = y;
								curNind[2] = z;

								if ( m_voxelID->GetPixel(curNind)>0 )
								{
									if(  m_labelingImage->GetPixel(curNind)>=0 )
									{
										adjlabel = m_labelingImage->GetPixel(curNind);
									}
								}
							}
						}

						xvec.push_back(w);
						yvec.push_back(h);
						zvec.push_back(d);
						int count = 1;
						for( int c=0; c<count; c++ )
						{
							for( int n=0; n<10; n++ )
							{
								int x = xvec[c] + dx10[n];
								int y = yvec[c] + dy10[n];
								int z = zvec[c] + dz10[n];

								if( (x>=0 && x<m_ImageSize[0]) && (y>=0 && y<m_ImageSize[1]) && (z>=0 && z<m_ImageSize[2]))
								{
									VectorIndexType curNind;
									curNind[0] = x;
									curNind[1] = y;
									curNind[2] = z;

									if ( m_voxelID->GetPixel(curNind)>0 )
									{
										if( 0>m_labelingImage->GetPixel(curNind) && m_OutputImage->GetPixel(curInd)==m_OutputImage->GetPixel(curNind) )
										{
											xvec.push_back(x);
											yvec.push_back(y);
											zvec.push_back(z);
											m_labelingImage->SetPixel( curNind,lab );
											
											count++;
										}
									}
								}
							}
						}
						//-------------------------------------------------------
						// If segment size is less than a limit, assign an
						// adjacent label found before, and decrement label count.
						//-------------------------------------------------------
						if(count <= m_minSize)
						{
							for( int c=0; c<count;c++ )
							{
								VectorIndexType ind;
								ind[0] = xvec[c];
								ind[1] = yvec[c];
								ind[2] = zvec[c];

								m_labelingImage->SetPixel( ind,adjlabel );
							}
							lab--;
						}

						//--------------------------------------------------------
						lab++;

						xvec.clear();
						yvec.clear();
						zvec.clear();
					}
				}
			}
		}
	}

	m_labNum = lab;
	//------------------
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::DoSupervoxelSegmentation()
{
	// extract features
	FeaExtraction();
	
	vector<float> nullVecFt;
	m_SeedCentroids.resize( m_feaNum+TDimension,nullVecFt );

	// get the initial centroids here
	GenerateSeeds();
	// perform SLIC
	PerformSLIC( m_compactness );
	// enforce label connectivity
	EnforceSupervoxelLabelConnectivity();
}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::GetSvxlFea( string imgName )
{
	if (m_logger_or_not)
		m_logger->Write("Get super-voxel feature ...");
	else
		std::cout << "Get super-voxel feature" << std::endl;

	typename ScalarImageReaderType::Pointer ImageReader = ScalarImageReaderType::New();
	ImageReader->SetFileName( imgName );
	ImageReader->Update();

	float minIntensity = FLT_MAX;
	float maxIntensity = 0;
	float minGraM = FLT_MAX;
	float maxGraM = 0;
	ScalarIteratorType intIter(ImageReader->GetOutput(), ImageReader->GetOutput()->GetRequestedRegion());
	ScalarIteratorType graIter(m_GraMagnitude, m_GraMagnitude->GetRequestedRegion());
	intIter.GoToBegin();
	graIter.GoToBegin();
	while(!intIter.IsAtEnd())
	{
		float curInt = intIter.Get();
		float curGra = graIter.Get();

		if (curInt>maxIntensity)
			maxIntensity = curInt;
		if (curInt<minIntensity)
			minIntensity = curInt;
		if (curGra>maxGraM)
			maxGraM = curGra;
		if (curGra<minGraM)
			minGraM = curGra;

		++intIter;
		++graIter;
	}

	typename LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  	labelStatisticsImageFilter->SetLabelInput( m_labelingImage );
  	labelStatisticsImageFilter->SetInput( ImageReader->GetOutput() );
  	labelStatisticsImageFilter->SetHistogramParameters( 32, minIntensity, maxIntensity );
  	labelStatisticsImageFilter->Update();

  	int numLab = labelStatisticsImageFilter->GetNumberOfLabels();

	//if (m_logger_or_not)
	//	m_logger->Write("	Number of labels: " + to_string(numLab));
	//else
 // 		cout << "	Number of labels: " << numLab << endl;
 	
	m_SvxFea.set_size(numLab-1, 10); // exclude label 0
	m_SvxFea.fill(0.0f);
    
  	typedef typename LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  	typedef typename LabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;
 	typename ValidLabelValuesType::const_iterator vIt;

  	for( vIt=labelStatisticsImageFilter->GetValidLabelValues().begin();
    	 vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
      	 ++vIt )
    {
    	if ( labelStatisticsImageFilter->HasLabel(*vIt) )
      	{
      		LabelPixelType labelValue = *vIt;

      		if ( labelValue>0 )
      		{
      			m_SvxFea(labelValue-1, 0) = labelStatisticsImageFilter->GetMinimum( labelValue );
      			m_SvxFea(labelValue-1, 1) = labelStatisticsImageFilter->GetMaximum( labelValue );
      			m_SvxFea(labelValue-1, 2) = labelStatisticsImageFilter->GetMedian( labelValue );
      			m_SvxFea(labelValue-1, 3) = labelStatisticsImageFilter->GetMean( labelValue );
      			m_SvxFea(labelValue-1, 4) = labelStatisticsImageFilter->GetSigma( labelValue );
      			//m_SvxFea(labelValue-1, 5) = labelStatisticsImageFilter->GetVariance( labelValue );
      		}
      	}
    }

    labelStatisticsImageFilter->SetInput( m_GraMagnitude );
  	labelStatisticsImageFilter->SetHistogramParameters( 32, minGraM, maxGraM );
  	labelStatisticsImageFilter->Update();

  	for( vIt=labelStatisticsImageFilter->GetValidLabelValues().begin();
    	 vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
      	 ++vIt )
    {
    	if ( labelStatisticsImageFilter->HasLabel(*vIt) )
      	{
      		LabelPixelType labelValue = *vIt;

      		if ( labelValue>0 )
      		{
      			m_SvxFea(labelValue-1, 5) = labelStatisticsImageFilter->GetMinimum( labelValue );
      			m_SvxFea(labelValue-1, 6) = labelStatisticsImageFilter->GetMaximum( labelValue );
      			m_SvxFea(labelValue-1, 7) = labelStatisticsImageFilter->GetMedian( labelValue );
      			m_SvxFea(labelValue-1, 8) = labelStatisticsImageFilter->GetMean( labelValue );
      			m_SvxFea(labelValue-1, 9) = labelStatisticsImageFilter->GetSigma( labelValue );
      			//m_SvxFea(labelValue-1, 10) = labelStatisticsImageFilter->GetVariance( labelValue );
      		}
      	}
    }

    m_SvxFea.normalize_columns();
    m_SvxFea *= 10.0; 
    //for ( unsigned int sri=0; sri<m_SvxFea.rows(); sri++ )
    //{
    //    cout << sri+1 << ": ";
    //  	for ( unsigned int sfi=0; sfi<m_SvxFea.columns(); sfi++)
	//      	cout << m_SvxFea(sri, sfi) << "  ";
	//    cout << endl;    	
    //}

}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::DoKmeans( int numK, string refImg )
{
	if (m_logger_or_not)
		m_logger->Write("Kmeans ...");
	else
		std::cout << "Kmeans ..." << std::endl;
	
	int numSample = m_SvxFea.rows();
	int numFea = m_SvxFea.columns();

	Mat dataMat(numSample, numFea, CV_32F);
	Mat labels;
	Mat centers;

	for ( unsigned int sri=0; sri<m_SvxFea.rows(); sri++ )
    {
      	for ( unsigned int sfi=0; sfi<m_SvxFea.columns(); sfi++)
        {
            dataMat.at<float>(sri,sfi) = m_SvxFea(sri,sfi);
            //cout << dataMat.at<float>(sri,sfi) << "  ";
        }
        //cout << endl;
    }	

	kmeans(dataMat, numK, labels, TermCriteria(TermCriteria::EPS+TermCriteria::COUNT, 1000, 0.0001), 5, KMEANS_PP_CENTERS, centers);
	//for ( int i=0; i<numSample; i++)
	//	cout << labels.at<int>(i) << endl;

	// output Kmeans results
	m_KmeanLabelingImage = ScalarLabelImageType::New();
	m_KmeanLabelingImage->SetOrigin( m_OutputImage->GetOrigin() );
	m_KmeanLabelingImage->SetSpacing( m_OutputImage->GetSpacing() );
	m_KmeanLabelingImage->SetRegions( m_OutputImage->GetRequestedRegion() );
	m_KmeanLabelingImage->SetDirection( m_OutputImage->GetDirection() );
	m_KmeanLabelingImage->Allocate(); 
	m_KmeanLabelingImage->FillBuffer( -1 );

	ScalarLabelImageIteratorType Siter( m_labelingImage,m_labelingImage->GetRequestedRegion() );
	ScalarLabelImageIteratorType Oiter( m_KmeanLabelingImage,m_KmeanLabelingImage->GetRequestedRegion() );

	Siter.GoToBegin();
	Oiter.GoToBegin();
	while(!Siter.IsAtEnd())
	{
		//if ( Siter.Get()==0 )
		//	Oiter.Set( 0 );
		//else
		//	Oiter.Set( labels.at<int>(Siter.Get()-1) + 1 );
		if ( Siter.Get()<=0 )
			Oiter.Set( 0 );
		else
			Oiter.Set( labels.at<int>(Siter.Get()-1) + 1 );

		++Siter;
		++Oiter;
	}

	// pick the region containing lung field
	typename ScalarImageReaderType::Pointer ImageReader = ScalarImageReaderType::New();
	ImageReader->SetFileName( refImg );
	ImageReader->Update();

	typename LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
  	labelStatisticsImageFilter->SetLabelInput( m_KmeanLabelingImage );
  	labelStatisticsImageFilter->SetInput( ImageReader->GetOutput() );
  	labelStatisticsImageFilter->Update();

  	typedef typename LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
  	typedef typename LabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;
 	typename ValidLabelValuesType::const_iterator vIt;

  	int numLab = labelStatisticsImageFilter->GetNumberOfLabels();

 	vnl_vector< float > km_img_val(numLab-1, 0.0f);
 	vnl_vector< int > km_img_idx(numLab-1, 0);

  	for( vIt=labelStatisticsImageFilter->GetValidLabelValues().begin();
    	 vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
      	 ++vIt )
    {
    	if ( labelStatisticsImageFilter->HasLabel(*vIt) )
      	{
      		LabelPixelType labelValue = *vIt;

      		if ( labelValue>0 )
      		{
      			km_img_val(labelValue-1) = labelStatisticsImageFilter->GetMean( labelValue );
      			km_img_idx(labelValue-1) = labelValue;
      		}
      	}
    }
    unsigned int min_ind = km_img_val.arg_min();

	Oiter.GoToBegin();
	while(!Oiter.IsAtEnd())
	{
		if ( Oiter.Get()==km_img_idx(min_ind) )
			Oiter.Set( 1 );
		else
			Oiter.Set( 0 );
		
		++Oiter;
	}

	//
	typename ConnectedComponentType::Pointer ConnectedComponentFilter = ConnectedComponentType::New();
	ConnectedComponentFilter->SetInput( m_KmeanLabelingImage );
	// Relabel the components in order of size
	typename RelabelType::Pointer Relabeler = RelabelType::New();
	Relabeler->SetInput( ConnectedComponentFilter->GetOutput() );
	Relabeler->Update();

	int numRoi = Relabeler->GetNumberOfObjects();
	float max_vol = Relabeler->GetSizeOfObjectsInPixels()[0];
	if (numRoi>1)
	{
		ScalarLabelImageIteratorType rliter(Relabeler->GetOutput(), Relabeler->GetOutput()->GetRequestedRegion());

		Oiter.GoToBegin();
		rliter.GoToBegin();
		while(!Oiter.IsAtEnd())
		{
			LabelPixelType cur_lab = rliter.Get();
			float cur_vol = Relabeler->GetSizeOfObjectsInPixels()[cur_lab-1];
			if (cur_vol/max_vol < 0.25)
				Oiter.Set(0);

			++Oiter;
			++rliter;
		}
	}

}

template < class TPixelType, unsigned int TDimension, unsigned int TInputImageNum >
void SBRT_LungField< TPixelType, TDimension, TInputImageNum >::WriteLabel(string outBaseName)// , int maskAv )
{
  if (m_logger_or_not)
    m_logger->Write("Write label");
  else
    std::cout << "Write label" << std::endl;

  // fg.nii.gz
  //if (maskAv==0)
  //{
  //	string fgImgName = outBaseName + "_fg.nii.gz";
    //
  //	typename ScalarImageWriterType::Pointer scalarWriter = ScalarImageWriterType::New();
  //	scalarWriter->SetFileName( fgImgName );
  //	scalarWriter->SetInput( m_Mask );
  //	scalarWriter->Update();
  //}


  typename ScalarLabelImageWriterType::Pointer labelWriter = ScalarLabelImageWriterType::New();
  //string slicImageName = outBaseName + "_slic.nii.gz";
  // slic.nii.gz
  //ScalarLabelImageIteratorType Siter( m_labelingImage,m_labelingImage->GetRequestedRegion() );
  //Siter.GoToBegin();
  //while(!Siter.IsAtEnd())
  //{
  //	if ( Siter.Get()<0 )
  //		Siter.Set( 0 );
  //	else
  //		Siter.Set( Siter.Get()+1 );
  //
  //	++Siter;
  //}

  //labelWriter->SetFileName( slicImageName );
  //labelWriter->SetInput( m_labelingImage );
  //labelWriter->Update();

  ScalarIteratorType Miter(m_Mask, m_Mask->GetRequestedRegion());
  ScalarLabelImageIteratorType Oiter(m_KmeanLabelingImage, m_KmeanLabelingImage->GetRequestedRegion());

  Miter.GoToBegin();
  Oiter.GoToBegin();
  while (!Miter.IsAtEnd())
  {
    if (Oiter.Get() > 0)
      Oiter.Set(2);
    else if (Miter.Get() > 0)
      Oiter.Set(3);

    ++Miter;
    ++Oiter;
  }

  string lfImageName = outBaseName + "_lf.nii.gz";
  // lf.nii.gz
  labelWriter->SetFileName(lfImageName);
  labelWriter->SetInput(m_KmeanLabelingImage);
  labelWriter->Update();
}

