#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkExtractImageFilter.h"
#include "itkCSVArray2DFileReader.h"

typedef itk::Image<float, 4> PerfusionImageType;
typedef itk::Image<float, 3> ImageType;
std::string inputImageFile, outputImageFile;
typedef itk::CSVArray2DFileReader<double> CSVFileReaderType;
typedef vnl_matrix<double> MatrixType;
using VariableLengthVectorType = itk::VariableLengthVector< double >; // TBD: move from itk::matrix to vnl::matrix for simplicity
using VariableSizeMatrixType = itk::VariableSizeMatrix< double >; // TBD: move from itk::matrix to vnl::matrix for simplicity

//interpolation function for one vector
void interp1(int *x, int x_tam, double *y, int *xx, int xx_tam, double *yy)
{
  double *dx, *dy, *slope, *intercept, *elementoMasProximo, *xD;
  int i, *indiceEnVector;

  dx = (double *)calloc(x_tam - 1, sizeof(double));
  dy = (double *)calloc(x_tam - 1, sizeof(double));
  slope = (double *)calloc(x_tam - 1, sizeof(double));
  intercept = (double *)calloc(x_tam - 1, sizeof(double));
  indiceEnVector = (int *)malloc(sizeof(int));
  elementoMasProximo = (double *)malloc(sizeof(double));
  xD = (double *)calloc(x_tam, sizeof(double));

  for (i = 0; i<x_tam; i++) {
    xD[i] = x[i];
  }

  for (i = 0; i < x_tam; i++) {
    if (i<x_tam - 1) {
      dx[i] = x[i + 1] - x[i];
      dy[i] = y[i + 1] - y[i];
      slope[i] = dy[i] / dx[i];
      intercept[i] = y[i] - x[i] * slope[i];
    }
    else {
      dx[i] = dx[i - 1];
      dy[i] = dy[i - 1];
      slope[i] = slope[i - 1];
      intercept[i] = intercept[i - 1];
    }
  }

 /* for (i = 0; i < xx_tam; i++) {
    encuentraValorMasProximo(xx[i], xD, x_tam, x_tam, elementoMasProximo, indiceEnVector);
    yy[i] = slope[*indiceEnVector] * xx[i] + intercept[*indiceEnVector];*/
  //}
}


template< class ImageType, class PerfusionImageType>
typename ImageType::Pointer GetOneImageVolume(typename PerfusionImageType::Pointer perfImagePointerNifti, int index)
{
  typename PerfusionImageType::RegionType region = perfImagePointerNifti->GetLargestPossibleRegion();
  typename PerfusionImageType::IndexType regionIndex;
  typename PerfusionImageType::SizeType regionSize;
  regionSize[0] = region.GetSize()[0];
  regionSize[1] = region.GetSize()[1];
  regionSize[2] = region.GetSize()[2];
  regionSize[3] = 0;
  regionIndex[0] = 0;
  regionIndex[1] = 0;
  regionIndex[2] = 0;
  regionIndex[3] = index;
  typename PerfusionImageType::RegionType desiredRegion(regionIndex, regionSize);

  typedef itk::ExtractImageFilter< PerfusionImageType, ImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(desiredRegion);
  filter->SetInput(perfImagePointerNifti);

  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();
  return filter->GetOutput();
}

int main(int argc, char** argv)
{
  cbica::CmdParser parser(argc, argv);

  //parser.addRequiredParameter("i", "inputImage", cbica::Parameter::FILE, "NIfTI", "Input Image for processing");
  //parser.addRequiredParameter("o", "outputImage", cbica::Parameter::FILE, "NIfTI", "Output Image for processing");

  //if (parser.isPresent("i"))
  //{
  //  parser.getParameterValue("i", inputImageFile);
  //}
  //if (parser.isPresent("o"))
  //{
  //  parser.getParameterValue("o", outputImageFile);
  //}

  //auto inputImageInfo = cbica::ImageInfo(inputImageFile);

  //switch (inputImageInfo.GetImageDimensions())
  //{
  //case 4:
  //{
    std::string inputFolder = "E:/SoftwareDevelopmentProjects/PerfusionAlignmentRelatedMaterial/PerfImages";
    std::string inputdicomFolder = "E:/SoftwareDevelopmentProjects/PerfusionAlignmentRelatedMaterial/Dicoms";
    std::string outputFolder = "E:/SoftwareDevelopmentProjects/PerfusionAlignmentRelatedMaterial";


    //read the corresponding perfusion images

    std::vector<std::string> subjectNames = cbica::filesInDirectory(inputFolder);
    std::sort(subjectNames.begin(), subjectNames.end());
    std::vector<std::vector<double>> BrainsAllPatients;
    subjectNames = cbica::filesInDirectory(inputFolder);
    std::sort(subjectNames.begin(), subjectNames.end());
    
    for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
    {
      std::cout << sid << ": " << subjectNames[sid] << std::endl;
      //reading the given perfusion image
      typename PerfusionImageType::Pointer perfImagePointerNifti;
      perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(subjectNames[sid]);

      //Get standard deviation 
      typename ImageType::Pointer StdDev = GetOneImageVolume<ImageType, PerfusionImageType>(perfImagePointerNifti, 0);
      for (unsigned int x = 0; x < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[0]; x++)
        for (unsigned int y = 0; y < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[1]; y++)
          for (unsigned int z = 0; z < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[2]; z++)
          {
            typename PerfusionImageType::IndexType index4D;
            index4D[0] = x;
            index4D[1] = y;
            index4D[2] = z;

            typename ImageType::IndexType index3D;
            index3D[0] = x;
            index3D[1] = y;
            index3D[2] = z;

            //---------------------------------------mean from 1-10------------------------------------
            std::vector<double> values;
            for (unsigned int k = 0; k < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3]; k++)
            {
              index4D[3] = k;
              values.push_back(perfImagePointerNifti->GetPixel(index4D));
            }
            //calcualte standard deviation of values
            double temp = 0.0;
            double stddev = 0;
            for (unsigned int sampleNo = 0; sampleNo < values.size(); sampleNo++)
              temp += values[sampleNo];
            double mean = temp / values.size();
            temp = 0.0;
            for (unsigned int sampleNo = 0; sampleNo < values.size(); sampleNo++)
              temp += (values[sampleNo] - mean)*(values[sampleNo] - mean);
            stddev = std::sqrt(temp / (values.size() - 1));
            StdDev.GetPointer()->SetPixel(index3D, stddev);
          }
      //calcualte average brain at all time points
      std::vector<double> averagebrain;
      for (int i = 0; i < 100; i++)
        averagebrain.push_back(0);
      for (unsigned int timePoints = 0; timePoints < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3]; timePoints++)
      {
        double summationvalues = 0;
        double counter = 0;
        for (unsigned int x = 0; x < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[0]; x++)
        {
          for (unsigned int y = 0; y < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[1]; y++)
          {
            for (unsigned int z = 0; z < perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[2]; z++)
            {
              typename PerfusionImageType::IndexType index4D;
              index4D[0] = x;
              index4D[1] = y;
              index4D[2] = z;
              index4D[3] = timePoints;

              typename ImageType::IndexType index3D;
              index3D[0] = x;
              index3D[1] = y;
              index3D[2] = z;

              if (StdDev->GetPixel(index3D) > 0)
              {
                summationvalues = summationvalues + perfImagePointerNifti->GetPixel(index4D);
                counter = counter + 1;
              }
            }
          }
        }
        averagebrain[timePoints] = summationvalues / counter;
      }
      BrainsAllPatients.push_back(averagebrain);
    }
    std::ofstream myfile;
    myfile.open(outputFolder + "/AverageBrains.csv");
    for (unsigned int index1 = 0; index1 < BrainsAllPatients.size(); index1++)
    {
      for (unsigned int index2 = 0; index2 < BrainsAllPatients[0].size(); index2++)
      {
        if (index2 == 0)
          myfile << std::to_string(BrainsAllPatients[index1][index2]);
        else
          myfile << "," <<std::to_string(BrainsAllPatients[index1][index2]);
      }
      myfile << "\n";
    }
    myfile.close();




   //read the tags from the dicom data
    std::vector<double> echoTime;
    subjectNames = cbica::filesInDirectory(inputdicomFolder);
    std::sort(subjectNames.begin(), subjectNames.end());
    for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
    {
      //reading dicom tags of the folders associated with these files
      using PixelType = signed short;
      constexpr unsigned int Dimension = 2;
      using ImageType = itk::Image< PixelType, Dimension >;
      using ReaderType = itk::ImageFileReader< ImageType >;
      ReaderType::Pointer reader = ReaderType::New();

      using ImageIOType = itk::GDCMImageIO;
      ImageIOType::Pointer dicomIO = ImageIOType::New();
      reader->SetFileName(subjectNames[sid]);
      reader->SetImageIO(dicomIO);
      reader->Update();

      using DictionaryType = itk::MetaDataDictionary;
      const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
      using MetaDataStringType = itk::MetaDataObject< std::string >;
      auto itr = dictionary.Begin();
      auto end = dictionary.End();
      std::string entryId = "0010|0010";
      auto tagItr = dictionary.Find(entryId);
      if (tagItr != end)
      {
        MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(tagItr->second.GetPointer());
        if (entryvalue)
        {
          //tagvalue::string tagvalue = entryvalue->GetMetaDataObjectValue();
          echoTime.push_back(0);
        }
      }

      std::cout << sid << ": " << subjectNames[sid] << std::endl;
    }
    std::ofstream myfile;
    myfile.open(outputFolder + "/EchoTime.csv");
    for (unsigned int index1 = 0; index1 < echoTime.size(); index1++)
    {
      myfile << std::to_string(echoTime[index1]) << ",";
      myfile << "\n";
    }
    myfile.close();



    //read csv files with the tags and with the defined cut off points
    std::vector<double> Tags;
    CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
    MatrixType dataMatrix;
    try
    {
      reader->SetFileName(outputFolder + "/EchoTime.csv");
      reader->SetFieldDelimiterCharacter(',');
      reader->HasColumnHeadersOff();
      reader->HasRowHeadersOff();
      reader->Parse();
      dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

      for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      {
        for (unsigned int j = 0; j < dataMatrix.cols(); j++)
          Tags.push_back(dataMatrix(i, j));
      }
    }
    catch (const std::exception& e1)
    {
      return false;
    }

    //read curve of the following patients
    std::vector<std::vector<double>> RevisedCurves;
    CSVFileReaderType::Pointer reader = CSVFileReaderType::New();
    MatrixType dataMatrix;
    try
    {
      reader->SetFileName(outputFolder + "/AverageBrains.csv");
      reader->SetFieldDelimiterCharacter(',');
      reader->HasColumnHeadersOff();
      reader->HasRowHeadersOff();
      reader->Parse();
      dataMatrix = reader->GetArray2DDataObject()->GetMatrix();

      for (unsigned int i = 0; i < dataMatrix.rows(); i++)
      {
        std::vector<double> innervector;
        for (unsigned int j = 0; j < dataMatrix.cols(); j++)
          innervector.push_back(dataMatrix(i, j));
        RevisedCurves.push_back(innervector);
      }
    }
    catch (const std::exception& e1)
    {
      return false;
    }

    //calculater different statistics of curves
    std::vector<std::vector<double>> NormalizedCurves;
    std::vector<double> min_curve;
    std::vector<double> max_curve;
    std::vector<double> drop_curve;
    std::vector<double> base_curve;
    for (int index = 0; index < RevisedCurves.size(); index++)
    {
      std::vector<double> currentcurve = RevisedCurves[index];
      min_curve.push_back(*std::min_element(std::begin(currentcurve), std::end(currentcurve)));
      max_curve.push_back(*std::max_element(std::begin(currentcurve), std::end(currentcurve)));

      for (int index2 = 0; index2 < currentcurve.size(); index2++)
        currentcurve[index2] = (currentcurve[index2] * 100) / (max_curve[index] - min_curve[index]);

       NormalizedCurves.push_back(currentcurve);
       
      for (int index2 = 0; index2 < currentcurve.size(); index2++)
        currentcurve[index2] = currentcurve[index2] - max_curve[index] + 300;
      
      std::vector<double>::iterator result = std::min_element(std::begin(currentcurve), std::end(currentcurve));
      drop_curve.push_back(std::distance(std::begin(currentcurve), result));
      max_curve.push_back(*std::max_element(std::begin(currentcurve), std::end(currentcurve)));
    }


    std::string inputFolder = "E:/SoftwareDevelopmentProjects/PerfusionAlignmentRelatedMaterial/PerfImages";
    for (unsigned int sid = 0; sid < subjectNames.size(); sid++)
    {
      typename PerfusionImageType::Pointer perfImagePointerNifti;
      perfImagePointerNifti = cbica::ReadImage<PerfusionImageType>(subjectNames[sid]);
      double end_value = std::round((echoTime[sid]*2/1000) * perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3]);
      std::vector<double> xq_array;
      std::vector<double> x_array;
      std::vector<double> v;
      for (int index = 0; index < end_value; index++)
        xq_array.push_back(index);
      for (int index = 0; index < end_value; index=index+ (echoTime[sid] * 2 / 1000))
        x_array.push_back(index);
      //make a check on the length of x-array if greater than v


      auto outputSize = perfImagePointerNifti->GetLargestPossibleRegion().GetSize();
      auto outputSpacing = perfImagePointerNifti->GetSpacing();
      outputSize[3] = outputSize[3] * outputSpacing[3];
      outputSpacing[3] = 1 / (echoTime[sid] * 2 / 1000);
      auto resampler = itk::ResampleImageFilter< PerfusionImageType, PerfusionImageType >::New();
      resampler->SetInput(perfImagePointerNifti);
      resampler->SetSize(outputSize);
      resampler->SetOutputSpacing(outputSpacing);
      resampler->SetOutputOrigin(perfImagePointerNifti->GetOrigin());
      resampler->SetOutputDirection(perfImagePointerNifti->GetDirection());
      resampler->SetOutputStartIndex(perfImagePointerNifti->GetLargestPossibleRegion().GetIndex());
      resampler->SetTransform(itk::IdentityTransform< double, PerfusionImageType::ImageDimension >::New());
      resampler->UpdateLargestPossibleRegion();

    }
 

    //break;

  return EXIT_SUCCESS;
}