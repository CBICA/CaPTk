#pragma once

#include <sstream>
#include "SBRT_Analysis.h"

template < class TPixelType, unsigned int TDimension >
SBRT_Analysis< TPixelType, TDimension >::SBRT_Analysis()
{
	m_roiLabel = 1;
	m_numRoi = 1;
    
    m_logger = NULL;
	m_logger_or_not = false;
}

template < class TPixelType, unsigned int TDimension >
SBRT_Analysis< TPixelType, TDimension >::~SBRT_Analysis()
{

}

template < class TPixelType, unsigned int TDimension >
int SBRT_Analysis< TPixelType, TDimension >::Initialization( std::string InputImageName, std::string InputLabelName )
{
	if (m_logger_or_not)
        m_logger->Write("Read images ...");
    else
        std::cout << "Read images ... " << std::endl;

	typename ScalarImageReaderType::Pointer InputReader = ScalarImageReaderType::New();
	InputReader->SetFileName( InputImageName );
	InputReader->Update();

	typename ScalarLabelImageReaderType::Pointer MaskReader = ScalarLabelImageReaderType::New();
	MaskReader->SetFileName( InputLabelName );
	MaskReader->Update();

	//m_inputImage = ScalarImageType::New();
	m_inputImage = InputReader->GetOutput();
	m_inputImage->DisconnectPipeline();

	//m_labelImage = ScalarImageType::New();
	m_labelImage = MaskReader->GetOutput();
	m_labelImage->DisconnectPipeline();

	for( unsigned int j=0;j<TDimension;j++ )
	{
		m_ImageVoxelSize.push_back( m_inputImage->GetSpacing()[j] );
	}

	m_ImageSize = m_inputImage->GetRequestedRegion().GetSize();

	//cout << "Read images successfully " << endl;
    
	//keep roi label only
	ScalarLabelImageIteratorType liter(m_labelImage, m_labelImage->GetRequestedRegion());
	liter.GoToBegin();
	while(!liter.IsAtEnd())
	{
		if (liter.Get()!=m_roiLabel)
			liter.Set(0);

		++liter;
	}

	//get connected component and relabel them according to their volume
	typename ConnectedComponentType::Pointer ConnectedComponentFilter = ConnectedComponentType::New();
	ConnectedComponentFilter->SetInput( m_labelImage );
	// Relabel the components in order of size
	typename RelabelType::Pointer Relabeler = RelabelType::New();
	Relabeler->SetInput( ConnectedComponentFilter->GetOutput() );
	Relabeler->Update();

	m_numRoi = Relabeler->GetNumberOfObjects();
	//cout << "  # of ROIs: " << m_numRoi << endl;

	m_relabImage = Relabeler->GetOutput();
	m_relabImage->DisconnectPipeline();

	m_relabScalarImage = ScalarImageType::New();
	m_relabScalarImage->SetOrigin( m_relabImage->GetOrigin() );
	m_relabScalarImage->SetSpacing( m_relabImage->GetSpacing() );
	m_relabScalarImage->SetRegions( m_relabImage->GetRequestedRegion() );
	m_relabScalarImage->SetDirection( m_relabImage->GetDirection() );
	m_relabScalarImage->Allocate(); 
	m_relabScalarImage->FillBuffer( -1 );

	ScalarLabelImageIteratorType rliter(m_relabImage, m_relabImage->GetRequestedRegion());
	ScalarIteratorType rlsiter(m_relabScalarImage, m_relabScalarImage->GetRequestedRegion());
	rliter.GoToBegin();
	rlsiter.GoToBegin();
	while(!rliter.IsAtEnd())
	{
		rlsiter.Set(rliter.Get());
		++rliter;
		++rlsiter;
	}

    return 0;
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::SetParameters( unsigned int roiLab )
{
	m_roiLabel = roiLab;
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::SetLogger( std::string logName )
{
    m_logger = new cbica::Logging;
	m_logger->UseNewFile(logName);
    m_logger_or_not = true;
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::FeaExtraction()
{
    if (m_logger_or_not)
        m_logger->Write("Extracting features ...");
    else
        std::cout << "Extracting features ..." << std::endl;

	CalcShapeFeature();
	CalcIntensityFeature();
	CalcGlcmFeature();
	CalcGlrlmFeature();

	m_fea_all.push_back(m_shape);
	m_fea_all.push_back(m_intensity);
	m_fea_all.push_back(m_glcm);
	m_fea_all.push_back(m_glrlm);
    
	m_feaName_all.push_back(m_shape_name);
	m_feaName_all.push_back(m_intensity_name);
	m_feaName_all.push_back(m_glcm_name);
	m_feaName_all.push_back(m_glrlm_name);
    
    for (unsigned int i=0; i<m_fea_all.size(); i++)
        for (unsigned int j=0; j<m_fea_all[i].size();j++)
        {
        	if (std::isnan(m_fea_all[i][j]))
        		m_fea_vec.push_back(0.0);
        	else
	            m_fea_vec.push_back(m_fea_all[i][j]);
        }
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::CalcShapeFeature()
{
	typename LabImg2LabMapType::Pointer i2l = LabImg2LabMapType::New();
	i2l->SetInput(m_relabImage);
	i2l->SetComputePerimeter(true);
	i2l->Update();

	LabelMapType *labelMap = i2l->GetOutput();
	//cout << "  # of label objects: " << labelMap->GetNumberOfLabelObjects() << endl;
	ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(0); // the first object
	//cout << "    label: " << labelObject->GetLabel() << endl;
	//ShapeLabelObjectType *labelObject = labelMap->GetLabelObject(1);

	m_shape_name.push_back("no_of_pix");
  	//m_shape_name.push_back("principleMomemt_1");
  	//m_shape_name.push_back("principleMomemt_2");
  	//m_shape_name.push_back("principleMomemt_3");
  	m_shape_name.push_back("Elongation");
  	m_shape_name.push_back("Perimeter");
  	m_shape_name.push_back("Roundness");
  	m_shape_name.push_back("Flatness");

  	m_shape.push_back(labelObject->GetNumberOfPixels());
  	//m_shape.push_back(labelObject->GetPrincipalMoments().operator[](0));
  	//m_shape.push_back(labelObject->GetPrincipalMoments().operator[](1));
  	//m_shape.push_back(labelObject->GetPrincipalMoments().operator[](2));
  	m_shape.push_back(labelObject->GetElongation());
  	m_shape.push_back(labelObject->GetPerimeter());
  	m_shape.push_back(labelObject->GetRoundness());
  	m_shape.push_back(labelObject->GetFlatness());
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::CalcIntensityFeature()
{
	typename LabImg2StatLabMapType::Pointer i2sl = LabImg2StatLabMapType::New();
	i2sl->SetInput1(m_relabImage);
	i2sl->SetInput2(m_inputImage);
	i2sl->Update();

	StatisticsLabelMapType *statLabMap = i2sl->GetOutput();

	//typename StatisticsRelabelMapType::Pointer StalRelabelMap = StatisticsRelabelMapType::New();
	//StalRelabelMap->SetInput(statLabMap);
	//StalRelabelMap->SetAttribute("Size");
	//StalRelabelMap->Update();

	//StatisticsLabelObjectType *stalObj = StalRelabelMap->GetOutput()->GetLabelObject(1);
	StatisticsLabelObjectType *stalObj = statLabMap->GetLabelObject( 1 );
	//cout << "    stat label: " << stalObj->GetLabel() << endl;
	
	m_intensity_name.push_back("Minimum");
	m_intensity_name.push_back("Maximum");
	m_intensity_name.push_back("Mean");
	m_intensity_name.push_back("StandardDeviation");
	m_intensity_name.push_back("Variance");
	m_intensity_name.push_back("Skewness");
	m_intensity_name.push_back("Kurtosis");

	m_intensity.push_back( stalObj->GetMinimum() );
	m_intensity.push_back( stalObj->GetMaximum() );
	m_intensity.push_back( stalObj->GetMean() );
	m_intensity.push_back( stalObj->GetStandardDeviation() );
	m_intensity.push_back( stalObj->GetVariance() );
	m_intensity.push_back( stalObj->GetSkewness() );
	m_intensity.push_back( stalObj->GetKurtosis() );
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::CalcGlcmFeature()
{
	unsigned int insideValue = 1;

	ScalarIteratorType iiter(m_inputImage,m_inputImage->GetRequestedRegion());
	ScalarLabelImageIteratorType miter(m_relabImage,m_relabImage->GetRequestedRegion());
	TPixelType intensityMin = 1e6;
	TPixelType intensityMax = -1e6;
	iiter.GoToBegin();
	miter.GoToBegin();
	while (!iiter.IsAtEnd())
	{
		if ( miter.Get()==insideValue )
		{
			if ( iiter.Get() < intensityMin )
				intensityMin = iiter.Get();
			if ( iiter.Get() > intensityMax )
				intensityMax = iiter.Get();
		}
		++iiter;
		++miter;
	}

	typename GlcmType::FeatureNameVectorPointer requestedFeatures = GlcmType::FeatureNameVector::New();
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::Energy);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::Entropy);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::InverseDifferenceMoment);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::Inertia);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::ClusterShade);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::ClusterProminence);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::Correlation);
  	requestedFeatures->push_back(GlcmType::TextureFeaturesFilterType::HaralickCorrelation);
  	
	typename GlcmType::Pointer Texture = GlcmType::New();
	Texture->SetInput( m_inputImage );
	Texture->SetMaskImage( m_relabScalarImage );
	Texture->SetInsidePixelValue( insideValue );
	Texture->SetPixelValueMinMax( intensityMin,intensityMax );
	Texture->SetNumberOfBinsPerAxis( 32 );
	Texture->SetRequestedFeatures(requestedFeatures);
	Texture->Update();

	m_glcm_name.push_back("GLCM_Energy");
	m_glcm_name.push_back("GLCM_Entropy");
	m_glcm_name.push_back("GLCM_InverseDifferenceMoment");
	m_glcm_name.push_back("GLCM_Inertia");
	m_glcm_name.push_back("GLCM_ClusterShade");
	m_glcm_name.push_back("GLCM_ClusterProminence");
	m_glcm_name.push_back("GLCM_Correlation");
	m_glcm_name.push_back("GLCM_HaralickCorrelation");

	const typename GlcmType::FeatureValueVector* tex = Texture->GetFeatureMeans();
	for( unsigned int i = 0; i < tex->size(); ++i )
	{
		m_glcm.push_back( (*tex)[i] );
	}
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::CalcGlrlmFeature()
{
	unsigned int insideValue = 1;

	ScalarIteratorType iiter(m_inputImage,m_inputImage->GetRequestedRegion());
	ScalarLabelImageIteratorType miter(m_relabImage,m_relabImage->GetRequestedRegion());
	TPixelType intensityMin = 1e6;
	TPixelType intensityMax = -1e6;
	iiter.GoToBegin();
	miter.GoToBegin();
	while (!iiter.IsAtEnd())
	{
		if ( miter.Get()==insideValue )
		{
			if ( iiter.Get() < intensityMin )
				intensityMin = iiter.Get();
			if ( iiter.Get() > intensityMax )
				intensityMax = iiter.Get();
		}
		++iiter;
		++miter;
	}

	typename GlrlmType::Pointer Texture = GlrlmType::New();
	Texture->SetInput( m_inputImage );
	Texture->SetMaskImage( m_relabScalarImage );
	Texture->SetInsidePixelValue( insideValue );
	Texture->SetPixelValueMinMax( intensityMin,intensityMax );
	Texture->SetNumberOfBinsPerAxis( 32 );
	Texture->Update();

	m_glrlm_name.push_back("GLRLM_ShortRunEmphasis");
	m_glrlm_name.push_back("GLRLM_LongRunEmphasis");
	m_glrlm_name.push_back("GLRLM_GreyLevelNonuniformity");
	m_glrlm_name.push_back("GLRLM_RunLengthNonuniformity");
	m_glrlm_name.push_back("GLRLM_LowGreyLevelRunEmphasis");
	m_glrlm_name.push_back("GLRLM_HighGreyLevelRunEmphasis");
	m_glrlm_name.push_back("GLRLM_ShortRunLowGreyLevelEmphasis");
	m_glrlm_name.push_back("GLRLM_ShortRunHighGreyLevelEmphasis");
	m_glrlm_name.push_back("GLRLM_LongRunLowGreyLevelEmphasis");
	m_glrlm_name.push_back("GLRLM_LongRunHighGreyLevelEmphasis");

	const typename GlrlmType::FeatureValueVector* tex = Texture->GetFeatureMeans();
	for( unsigned int i = 0; i < tex->size(); ++i )
	{
		m_glrlm.push_back( (*tex)[i] );
	}
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::OutputFeature( std::string oname )
{
    if (m_logger_or_not)
        m_logger->Write("Write the features into file ...");
    else
        std::cout << "Write the features into file ..." << std::endl;

	std::string ofea_file_bin = oname + ".bin";
	std::ofstream soutfile( ofea_file_bin, std::ios::app|std::ios::binary );
	
	std::string ofea_file_asc = oname + "_asc.txt";
	std::ofstream soutfileASC( ofea_file_asc, std::ios::app );

	std::string oname_file_asc = oname + "_feaName.txt";
	std::ofstream soutfileInfo( oname_file_asc, std::ios::app );
	
	unsigned int totalFeaSize = 0;
	for ( unsigned int i=0; i<m_fea_all.size(); i++ )
	{
		for ( unsigned int j=0; j<m_fea_all[i].size(); j++ )
		{
			float curValue = m_fea_all[i][j];
			soutfile.write( (char *)(&curValue), sizeof(float) );
			soutfileASC << curValue << std::endl;

			soutfileInfo << m_feaName_all[i][j] << std::endl;
		}

		totalFeaSize += m_fea_all[i].size();
	}	
		
	soutfile.close();
	soutfileASC.close();
	soutfileInfo.close();
}

template < class TPixelType, unsigned int TDimension >
void SBRT_Analysis< TPixelType, TDimension >::GetPredictedRisk( std::string metaName, std::string riskProjName )
{
    if (m_logger_or_not)
        m_logger->Write("Get predicted risks ...");
    else
        std::cout << "Get predicted risks ..." << std::endl;

    std::vector<std::vector<float> > meta_proj;
    std::string line;
    
    std::ifstream ifile(metaName, std::ios::in);
    while( getline(ifile, line) )
    {
        std::istringstream reader(line);
        std::vector<float> lineData;
        
        float val = 0;
        while( reader>>val )
            lineData.push_back(val);
        
        meta_proj.push_back(lineData);
    }
    //for (unsigned int i=0; i<meta_proj.size(); i++)
    //{
    //    for (unsigned int j=0; j<meta_proj[i].size(); j++)
    //        cout << meta_proj[i][j] << "\t";
    //    cout << endl;
    //}
    
    std::vector<std::vector<float> > risk_proj;
    std::ifstream rpfile(riskProjName, std::ios::in);
    while( getline(rpfile, line) )
    {
        std::istringstream rpreader(line);
        std::vector<float> rplineData;
        
        float rval = 0;
        while( rpreader>>rval )
            rplineData.push_back(rval);
        
        risk_proj.push_back(rplineData);
    }
    //for (unsigned int i=0; i<risk_proj.size(); i++)
    //{
    //    for (unsigned int j=0; j<risk_proj[i].size(); j++)
    //        cout << risk_proj[i][j] << "\t";
    //    cout << endl;
    //}
    
    // data normalization
    std::vector<float> norm_fea_vec(m_fea_vec.size(), 0);
    for (unsigned int i=0; i<m_fea_vec.size(); i++)
    {
        norm_fea_vec[i] = (m_fea_vec[i]-meta_proj[0][i]) / (meta_proj[1][i]-meta_proj[0][i]+1e-9);
    }
    
    // meta-feature
    std::vector<float> meta_fea(meta_proj.size()-2, 0);
    for (unsigned int i=0; i<meta_proj.size()-2; i++)
    {
        float mf_val = 0;
        for (unsigned int j=0; j<meta_proj[i+2].size(); j++)
            mf_val += meta_proj[i+2][j] * norm_fea_vec[j];
        meta_fea[i] = mf_val;
    }
    
    // prediction
    std::vector<std::string> tmpRiskName;
    tmpRiskName.push_back("survival");
    tmpRiskName.push_back("nodal failure");

    int numRisk = risk_proj.size() / 2;
    m_predRisk.resize(numRisk);
    for (int i = 0; i < numRisk; i++)
    {
        float tmp_risk = 0;
        for (unsigned int j=0; j<risk_proj[2*i].size(); j++)
            tmp_risk += risk_proj[2*i+1][j] * (meta_fea[j]-risk_proj[2*i][j]);
        m_predRisk[i] = exp(tmp_risk);
        
        if (m_logger_or_not)
            m_logger->Write("Predicted risk (" + tmpRiskName[i] + "): " + std::to_string( m_predRisk[i]) );
        else
            std::cout << "Predicted risk (" + tmpRiskName[i] + "): " + std::to_string( m_predRisk[i]) << std::endl;
    }
    //cout << endl;
}

template < class TPixelType, unsigned int TDimension >
float SBRT_Analysis< TPixelType, TDimension >::GetSurivalRisk()
{
  return m_predRisk[0];
}

template < class TPixelType, unsigned int TDimension >
float SBRT_Analysis< TPixelType, TDimension >::GetNodalFailureRisk()
{
  return m_predRisk[1];
}
