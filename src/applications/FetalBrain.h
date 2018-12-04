#pragma once

#include <iostream>
#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkOrientImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkLaplacianImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkLabelToRGBImageFilter.h"
#include "QuickView.h"
#include "itkBilateralImageFilter.h"
#include <itkLabelImageToShapeLabelMapFilter.h>
#include <itkShapeLabelObject.h>
#include <itkLabelMapToBinaryImageFilter.h>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMedianImageFilter.h"

#include <opencv2/opencv.hpp>
#include "QuickView.h"
#include "itkOpenCVImageBridge.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include "cbicaUtilities.h"
#include "cbicaITKSafeImageIO.h"
#include <QDir>

//#include "CAPTk.h"

typedef itk::ImageRegionIteratorWithIndex< ImageTypeFloat3D > IteratorType;

class Fetalbrain{

public:

  Fetalbrain()
  {

  }
  ~Fetalbrain()
  {
  }


  ImageTypeFloat3D::Pointer  apply_slicemask(ImageTypeFloat3D::Pointer fixedimagepointer, ImageTypeFloat3D::Pointer mask2d, std::string subject)
  {
    //Copy information if you need to keep the original image unharmed
    //Alternatively the filter may be run on the original image directly
    auto brainmask3d = ImageTypeFloat3D::New();;
    brainmask3d->CopyInformation(fixedimagepointer);
    brainmask3d->SetRequestedRegion(fixedimagepointer->GetRequestedRegion());
    brainmask3d->SetBufferedRegion(fixedimagepointer->GetBufferedRegion());
    brainmask3d->SetOrigin(fixedimagepointer->GetOrigin());
    brainmask3d->SetDirection(fixedimagepointer->GetDirection());
    brainmask3d->SetSpacing(fixedimagepointer->GetSpacing());
    brainmask3d->Allocate();
    brainmask3d->FillBuffer(0);


    ImageTypeFloat3D::IndexType start_index, end_index;
    std::vector< ImageTypeFloat3D::IndexType> index_vec;

    IteratorType imgiter(fixedimagepointer, fixedimagepointer->GetLargestPossibleRegion());
    IteratorType iter(brainmask3d, brainmask3d->GetLargestPossibleRegion());
    IteratorType it(mask2d, mask2d->GetLargestPossibleRegion());
    imgiter.GoToBegin();
    for (it.GoToBegin(), iter.GoToBegin(); !it.IsAtEnd(); ++it, ++imgiter)
    {
      if (it.Get() != 0)
      {
        index_vec.push_back(it.GetIndex());
        iter.Set(imgiter.Get());
      }
      else
      {
        iter.Set(0);
      }
    }
    if (index_vec.size() != 0)
    {
      start_index = index_vec[0];
      end_index = index_vec[index_vec.size() - 1];

    }

    if (start_index[2] == end_index[2])
    {
      int slicestart = start_index[2] - 3;
      int sliceend = start_index[2] + 3;
      for (int j = slicestart; j < sliceend; j++)
      {
        for (int i = 0; i < index_vec.size(); i++)
        {
          auto index = index_vec[i];
          index[2] = j;
          iter.SetIndex(index);
          imgiter.SetIndex(index);
          iter.Set(imgiter.Get());
        }
      }
    }
   // cbica::WriteImage(brainmask3d, subject + "_3d_slice.nii.gz");

    return brainmask3d;

  }

  double linear_features(ImageTypeFloat3D::Pointer maskimagepointer, std::map<std::string, float>& featurevec, double FOHR = 0, double d = 0, double slice_no = 0)
  {
    std::multimap< int, ImageTypeFloat3D::IndexType> index_vec;
    typedef itk::ImageRegionIteratorWithIndex< ImageTypeFloat3D > IteratorType;
    IteratorType it(maskimagepointer, maskimagepointer->GetLargestPossibleRegion());

    double a = 0, b = 0, c = 0;

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
      //index_vec.insert(std::pair<int, ImageTypeFloat3D::IndexType>(statit_cast<int>(it.Get()), it.GetIndex()));
      if (it.Get() == 2.0)
      {
        index_vec.insert(std::pair<int, ImageTypeFloat3D::IndexType>(1, it.GetIndex()));
      }
      else if (it.Get() == 3.0)
      {
        index_vec.insert(std::pair<int, ImageTypeFloat3D::IndexType>(2, it.GetIndex()));
      }
      else if (it.Get() == 4.0)
      {
        index_vec.insert(std::pair<int, ImageTypeFloat3D::IndexType>(3, it.GetIndex()));
      }
      else if (it.Get() == 5.0)
      {
        index_vec.insert(std::pair<int, ImageTypeFloat3D::IndexType>(4, it.GetIndex()));
      }
    }
    decltype(index_vec.equal_range(1)) range;
    for (auto i = index_vec.begin(); i != index_vec.end(); i = range.second)
    {
      // Get the range of the current key
      range = index_vec.equal_range(i->first);

      itk::Point<double, 2> p1;
      p1[0] = range.first->second[0];
      p1[1] = range.first->second[1];
      slice_no = range.first->second[2];

      auto pp = range.second;
      --pp;

      itk::Point<double, 2> p2;
      p2[0] = pp->second[0];
      p2[1] = pp->second[1];

      if (i->first == 1)
        a = p2.EuclideanDistanceTo(p1);
      if (i->first == 2)
        b = p2.EuclideanDistanceTo(p1);
      if (i->first == 3)
        c = p2.EuclideanDistanceTo(p1);
      if (i->first == 4)
        d = p2.EuclideanDistanceTo(p1);
    }
    FOHR = (a + c) / (2 * b);

    //write to file the linear features
    featurevec["FOHR"] = FOHR;
    featurevec["Mean_AD"] = d;

    return slice_no;
  }


  ImageTypeFloat3D::Pointer segment(ImageTypeFloat3D::Pointer fixedimagepointer, int sliceno, std::string subjectUnderConsideration)
  {
    auto writer = itk::ImageFileWriter< ImageTypeFloat3D >::New();
    //  normalize input image 
    typedef itk::NormalizeImageFilter< ImageTypeFloat3D, ImageTypeFloat3D >
      NormalizeFilterType;
    NormalizeFilterType::Pointer normalizeFilter = NormalizeFilterType::New();
    normalizeFilter->SetInput(fixedimagepointer);
    normalizeFilter->GetOutput()->SetRequestedRegion(fixedimagepointer->GetRequestedRegion());
    normalizeFilter->Update();
    ImageTypeFloat3D::Pointer normalized_img = normalizeFilter->GetOutput();
    //writer->SetInput(normalizeFilter->GetOutput());
    //writer->SetFileName(subjectUnderConsideration + "/normalize.nii.gz");
    //writer->Update();

    //********************************///
    ///smoothing of input image /////
    typedef itk::DiscreteGaussianImageFilter<ImageTypeFloat3D, ImageTypeFloat3D >  filterType;
    filterType::Pointer gaussianFilter = filterType::New();
    gaussianFilter->SetInput(normalizeFilter->GetOutput());
    gaussianFilter->SetVariance(1);
    gaussianFilter->Update();
    //writer->SetInput(gaussianFilter->GetOutput());
    //writer->SetFileName("E:\\fetal_brain\\testdata1\\149\\smooth.nii");
    //writer->Update();

    typedef itk::MedianImageFilter< ImageTypeFloat3D, ImageTypeFloat3D >edgeFilterType;
    edgeFilterType::Pointer medianFilter = edgeFilterType::New();
    medianFilter->SetInput(normalizeFilter->GetOutput());
    edgeFilterType::InputSizeType radius;
    radius.Fill(3);
    medianFilter->Update();

    //writer->SetInput(medianFilter->GetOutput());
    //writer->SetFileName(subjectUnderConsideration + "/smoothedge.nii.gz");
    //writer->Update();

    typedef itk::ThresholdImageFilter <ImageTypeFloat3D>
      ThresholdImageFilterType;

    ThresholdImageFilterType::Pointer thresholdFilter
      = ThresholdImageFilterType::New();
    thresholdFilter->SetInput(normalizeFilter->GetOutput());
    thresholdFilter->ThresholdBelow(11);
    thresholdFilter->Update();

    //writer->SetInput(thresholdFilter->GetOutput());
    //writer->SetFileName(subjectUnderConsideration + "/seg.nii");
    //writer->Update();





    typedef itk::Image< unsigned short, 3 > OutputImageType;
    typedef itk::ConnectedComponentImageFilter <ImageTypeFloat3D, OutputImageType >
      ConnectedComponentImageFilterType;

    ConnectedComponentImageFilterType::Pointer connected =
      ConnectedComponentImageFilterType::New();
    connected->SetInput(thresholdFilter->GetOutput());
    connected->Update();

    typedef itk::RGBPixel<unsigned char>         RGBPixelType;
    typedef itk::Image<RGBPixelType, 3>  RGBImageType;
    typedef itk::LabelToRGBImageFilter<OutputImageType, RGBImageType> RGBFilterType;

    RGBFilterType::Pointer rgbFilter = RGBFilterType::New();
    rgbFilter->SetInput(connected->GetOutput());
    rgbFilter->Update();
    typedef itk::ImageFileWriter<RGBImageType>RGBWriterType;
    RGBImageType::Pointer rgbtimage = rgbFilter->GetOutput();
    //writerrgb->SetInput(rgbFilter->GetOutput());
    //writerrgb->SetFileName(subjectUnderConsideration + "/compo.nii");
    //writerrgb->Update();

    typedef itk::LabelShapeKeepNObjectsImageFilter< OutputImageType > LabelShapeKeepNObjectsImageFilterType;
    LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
    labelShapeKeepNObjectsImageFilter->SetInput(connected->GetOutput());
    labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
    labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
    labelShapeKeepNObjectsImageFilter->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);

    typedef itk::RescaleIntensityImageFilter< OutputImageType, ImageTypeFloat3D > RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(itk::NumericTraits<float>::max());
    rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
    rescaleFilter->Update();
    auto ventimage = rescaleFilter->GetOutput();


    typedef itk::ImageRegionIteratorWithIndex<ImageTypeFloat3D> IteratorType;
      typedef itk::ImageRegionIteratorWithIndex<RGBImageType> RGBIteratorType;
    RGBIteratorType it_seg(rgbFilter->GetOutput(), rgbFilter->GetOutput()->GetLargestPossibleRegion());

    typedef itk::ImageRegionIteratorWithIndex<ImageTypeFloat3D> IteratorType1;
    IteratorType1 iter_ventri(ventimage, ventimage->GetLargestPossibleRegion());
    iter_ventri.GoToBegin();
    while (!iter_ventri.IsAtEnd())
    {
      it_seg.SetIndex(iter_ventri.GetIndex());
     
      RGBPixelType pixel = it_seg.Get();
      int pixvalue;
      if (!pixel.GetBlue() && !pixel.GetGreen() && !pixel.GetRed())
        pixvalue = 0;
      else if (pixel.GetBlue() || pixel.GetGreen() || pixel.GetRed())
        pixvalue = 1;
      if (pixvalue != 0 && iter_ventri.Get() != 0)
      {
       iter_ventri.Set(1);
      }
      if (pixvalue == 1 && iter_ventri.Get() == 0)
      {
        iter_ventri.Set(2);
      }

      ++iter_ventri;
    }


    writer->SetInput(ventimage);
    writer->SetFileName(subjectUnderConsideration + "/ventri.nii");
    writer->Update();
    return ventimage;
 }



 int calculatelabel(std::string subjectid, std::string gt_file)
 {


   std::map<std::string, float> grounthtruth;
   if (gt_file.empty())
   {
     std::cout << "Provide a groundtruth for training";
     return  0;
   }
   else
   {
     std::ifstream file(gt_file);
     std::string key, value;
     while (file.good())
     {
       getline(file, key, '\n');
       while (getline(file, key, ','))
       {
         getline(file, value);
         grounthtruth[key] = stof(value);
       }

     }
   }

   int label = grounthtruth.find(subjectid)->second;
   return label;
 }


 cv::Mat testsubject(std::map<std::string, float> &featurevec)
 {
   cv::Mat testfeatures(1, featurevec.size(), CV_32F);
   int j = 0;
   for (auto const& x : featurevec)
   {
     if (j < featurevec.size())
     {
       testfeatures.at<float>(0, j) = x.second;
       if (testfeatures.at<float>(0, j) != testfeatures.at<float>(0, j) || std::isinf(testfeatures.at<float>(0, j)))
         testfeatures.at<float>(0, j) = 0;
       j++;
     }
   }
   std::string trainedmodel;
   if (QDir(QApplication::applicationDirPath() + "/../data/").exists()) // packaged binary
   {
     trainedmodel = QApplication::applicationDirPath().toStdString() + "/../data/FetalModel.xml";
   }
   else if (QDir(QApplication::applicationDirPath() + "/../../data/").exists()) // running from project location
   {
     trainedmodel = QApplication::applicationDirPath().toStdString() + "/../../data/FetalModel.xml";
   }
   auto svmmodel = cv::Algorithm::load<cv::ml::SVM>(trainedmodel);
   cv::Mat predicted(1, 1, CV_32F);
   cv::normalize(testfeatures, testfeatures, cv::NORM_MINMAX, CV_32F);
   float value = svmmodel->predict(testfeatures, predicted, cv::ml::StatModel::RAW_OUTPUT);
   return predicted;

 }

 void training(std::string datadirectory, std::string outputdirectory,std::map<std::string, float>& featurevec)
 {
   std::vector< std::string > detectedSubjects = cbica::subdirectoriesInDirectory(datadirectory);

   std::string groundtruthfile = datadirectory + "/gt.csv";

   std::map<std::string, float> grounthtruth;
   if (groundtruthfile.empty())
   {
     std::cout << "Provide a groundtruth for training";
     return  ;
   }
   else
   {
     std::ifstream file(groundtruthfile);
     std::string key, value;
     while (file.good())
     {
       getline(file, key, '\n');
       while (getline(file, key, ','))
       {
         getline(file, value);
         grounthtruth[key] = stof(value);
       }

     }
   }
   std::vector<std::string> subj;
   for (size_t subject = 0; subject < detectedSubjects.size(); subject++)
   {
     if (grounthtruth.find(detectedSubjects[subject]) != grounthtruth.end())
     {
       subj.push_back(detectedSubjects[subject]);
     }

   }
   cv::Mat trainingfeatures(subj.size(), 18, CV_32F);
   cv::Mat traininglabels(subj.size(), 1, CV_32S);
   for (size_t subject = 0; subject < subj.size(); subject++)
   {
     const std::string subjectUnderConsideration = datadirectory + "/" + subj[subject] + "/";

       std::vector< std::string > filesPerSubject = cbica::filesInDirectory(subjectUnderConsideration);

       auto inputimage = ImageTypeFloat3D::New();
       auto mask = ImageTypeFloat3D::New();

       for (size_t j = 0; j < filesPerSubject.size(); j++)
       {
         FileNameParts  fileUnderConsideration = FileNameParts(filesPerSubject[j]);


         if (!cbica::fileExists(fileUnderConsideration.fullFileName))
         {
           std::cout << "Aborting:   Could not read " << fileUnderConsideration.fullFileName << std::endl;

         }
         if (fileUnderConsideration.base == detectedSubjects[subject] + "_origScan")
         {
           if (!cbica::fileExists(fileUnderConsideration.fullFileName))
           {
             std::cout << "Aborting:   Could not read " << fileUnderConsideration.fullFileName << std::endl;

           }
           else
             inputimage = cbica::ReadImage<ImageTypeFloat3D>(fileUnderConsideration.fullFileName);
         }
         if (fileUnderConsideration.base == detectedSubjects[subject] + "_mask")
         {
           if (!cbica::fileExists(fileUnderConsideration.fullFileName))
           {
             std::cout << "Aborting:   Could not read " << fileUnderConsideration.fullFileName << std::endl;

           }
           else
             mask = cbica::ReadImage<ImageTypeFloat3D>(fileUnderConsideration.fullFileName);
         }

       }

       //   std::cout << " Input Images read.." << std::endl;

       double  d = 0, FOHR = 0, slice_no = 0;
       //  ImageTypeFloat3D::Pointer processed_input = apply_slicemask(inputimage, slice_mask, "");
       slice_no = linear_features(mask, featurevec);

       // segment(processed_input, slice_no,  subjectUnderConsideration);
       Calculatefeatures(featurevec, slice_no, subjectUnderConsideration);


       int j = 0;

       for (auto const& x : featurevec)
       {
         if (j < featurevec.size())
         {
           trainingfeatures.at<float>(subject, j) = x.second;
           if (trainingfeatures.at<float>(subject, j) != trainingfeatures.at<float>(subject, j) || std::isinf(trainingfeatures.at<float>(subject, j)))
             trainingfeatures.at<float>(subject, j) = 0;
           j++;
         }
         if (j == featurevec.size())
         {
       
           traininglabels.at<int>(subject, 0) = grounthtruth.find(subj[subject])->second;
         }

       }

   
   
   }
     //normalize
     cv::normalize(trainingfeatures, trainingfeatures, cv::NORM_MINMAX, CV_32F);
     std::string filename1 = outputdirectory + "/training_features_norm.csv";
     std::ofstream myfile;
     myfile.open(filename1.c_str());
     myfile << cv::format(trainingfeatures, cv::Formatter::FMT_CSV) << std::endl;
     myfile.close();

     auto svm = cv::ml::SVM::create();
     svm->setKernel(cv::ml::SVM::LINEAR);
     svm->setType(cv::ml::SVM::C_SVC);
     svm->setC(10);
     //   svm->setGamma(0.5);
     cv::Mat weights = (cv::Mat_<double>(2, 1) << 1, 1.5);
     svm->setClassWeights(weights);
     svm->train(trainingfeatures, cv::ml::ROW_SAMPLE, traininglabels);

     std::string svmfilename = outputdirectory + "/FetalModel.xml";
     svm->save(svmfilename);
 
 }

 void Calculatefeatures(std::map<std::string, float>& featurevec, int sliceno, std::string subjectUnderConsideration)
 {
   ImageTypeFloat3D::Pointer ventimage;
   if (cbica::fileExists(subjectUnderConsideration + "/corrected_ventri.nii"))
   {
     ventimage = cbica::ReadImage<ImageTypeFloat3D>(subjectUnderConsideration + "/ventri.nii");
   }
   else
   {
     ventimage = cbica::ReadImage<ImageTypeFloat3D>(subjectUnderConsideration + "/ventri.nii");
     }

   ShapeFeatures(ventimage, sliceno, featurevec, "Vent");
   ShapeFeatures(ventimage, sliceno, featurevec, "Sas");
 }


 void ShapeFeatures(ImageTypeFloat3D::Pointer mask, int sliceno, std::map<std::string, float>& featurevec, std::string maskname)
 {
   double sum_pix = 0;
   typedef itk::ImageRegionIteratorWithIndex<ImageTypeFloat3D> IteratorType;
   IteratorType it_mask(mask, mask->GetLargestPossibleRegion());
   if (maskname == "Sas")
   {
     while (!it_mask.IsAtEnd())
     {
       if (it_mask.Get() == 2)
       {
         ++sum_pix;
       }
       ++it_mask;
     }
   }
   if (maskname == "Vent")
   {
     while (!it_mask.IsAtEnd())
     {
       if (it_mask.Get() == 1)
       {
         ++sum_pix;
       }
       ++it_mask;
     }
   }


   typedef  short LabelType;
   // 3d features
   typedef itk::Image< LabelType, 3 > OutputImageType;
   typedef itk::ShapeLabelObject< LabelType, 3 > ShapeLabelObjectType;
   typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

   typedef itk::ConnectedComponentImageFilter < ImageTypeFloat3D, OutputImageType > ConnectedComponentImageFilterType;
   typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType > I2LType;
   ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
   connected->SetInput(mask);
   connected->Update();


   I2LType::Pointer i2l = I2LType::New();
   i2l->SetInput(connected->GetOutput());
   i2l->SetComputePerimeter(true);
   i2l->Update();
   LabelMapType::Pointer labelMap = i2l->GetOutput();

   ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(0);
   if (maskname == "Vent"){
     featurevec["VentVolume"] = sum_pix;
     featurevec["MajorAxis"] = labelObject->GetPrincipalMoments().operator[](0);
     featurevec["MinorAxis"] = labelObject->GetPrincipalMoments().operator[](1);
     featurevec["VentPerimeter"] = labelObject->GetPerimeterOnBorder();
     featurevec["VentPerimeterRatio"] = labelObject->GetPerimeterOnBorderRatio() / sum_pix;
     featurevec["VentThickness"] = labelObject->GetFlatness();
   }
   if (maskname == "Sas")
   {
     featurevec["SasVolume"] = sum_pix;
     featurevec["SasPerimeterRatio"] = labelObject->GetPerimeterOnBorderRatio() / sum_pix;

   }
   featurevec["RatioVenttoSaS"] = featurevec["VentVolume"] / featurevec["SasVolume"];


   //calculating 2D features
   using ImageTypeFloat2D = itk::Image < float, 2 >;
   typedef itk::ExtractImageFilter<ImageTypeFloat3D, ImageTypeFloat2D> ExtractSlice;
   ExtractSlice::Pointer extractslice = ExtractSlice::New();
   extractslice->SetInput(mask);


   auto inputRegion = mask->GetLargestPossibleRegion();
   auto size = inputRegion.GetSize();
   size[2] = 0;

   auto start = inputRegion.GetIndex();
   start[2] = sliceno;

   ImageTypeFloat3D::RegionType desiredRegion;
   desiredRegion.SetSize(size);
   desiredRegion.SetIndex(start);
   extractslice->SetExtractionRegion(desiredRegion);
   extractslice->SetDirectionCollapseToSubmatrix();
   extractslice->Update();

   // 2d features
   typedef itk::Image< LabelType, 2 > OutputImageType2D;
   typedef itk::ShapeLabelObject< LabelType, 2 > ShapeLabelObjectType2D;
   typedef itk::LabelMap< ShapeLabelObjectType2D > LabelMapType2D;

   typedef itk::ConnectedComponentImageFilter < ImageTypeFloat2D, OutputImageType2D > ConnectedComponentImageFilterType2D;
   typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType2D, LabelMapType2D > I2LType2D;
   ConnectedComponentImageFilterType2D::Pointer connected2D = ConnectedComponentImageFilterType2D::New();
   connected2D->SetInput(extractslice->GetOutput());
   connected2D->Update();


   I2LType2D::Pointer i2l2D = I2LType2D::New();
   i2l2D->SetInput(connected2D->GetOutput());
   i2l2D->SetComputePerimeter(true);
   i2l2D->Update();
   LabelMapType2D::Pointer labelMap2D = i2l2D->GetOutput();

   ShapeLabelObjectType2D::Pointer labelObject2D = labelMap2D->GetNthLabelObject(0);
   typedef itk::ImageRegionIteratorWithIndex<ImageTypeFloat2D> IteratorType2D;
   IteratorType2D it_mask2D(extractslice->GetOutput(), extractslice->GetOutput()->GetLargestPossibleRegion());
   sum_pix = 0;
   if (maskname == "Vent"){
     while (!it_mask2D.IsAtEnd())
     {
       if (it_mask2D.Get() == 1)
       {
         ++sum_pix;
       }
       ++it_mask2D;
     }
   }
   if (maskname == "Sas"){
     while (!it_mask2D.IsAtEnd())
     {
       if (it_mask2D.Get() == 2)
       {
         ++sum_pix;
       }
       ++it_mask2D;
     }
   }


   if (maskname == "Vent"){
     featurevec["VentArea"] = sum_pix;
     // featurevec["MajorAxis"] = labelObject->GetPrincipalMoments().operator[](0);
     //featurevec["MinorAxis"] = labelObject->GetPrincipalMoments().operator[](1);
     featurevec["VentPerimeter2D"] = labelObject2D->GetPerimeterOnBorder();
     featurevec["VentPerimeterRatio2D"] = labelObject2D->GetPerimeterOnBorderRatio() / sum_pix;
     featurevec["VentThickness2D"] = labelObject2D->GetFlatness();
   }
   if (maskname == "Sas")
   {
     featurevec["SasArea"] = sum_pix;
     featurevec["SasPerimeterRatio2D"] = labelObject2D->GetPerimeterOnBorderRatio() / sum_pix;

   }
   featurevec["RatioVenttoSaS2D"] = featurevec["VentArea"] / featurevec["SasArea"];

 }


  
 void convertITKtoOpenCV( ImageTypeFloat3D::Pointer inputimage, cv::Mat &output)
 {
   auto region = inputimage->GetLargestPossibleRegion();
   auto size = region.GetSize();
   unsigned int w = static_cast< unsigned int >(size[0]);
   unsigned int h = static_cast< unsigned int >(size[1]);
   IteratorType iter(inputimage, inputimage->GetLargestPossibleRegion());
   iter.GoToBegin();
   output.create(cv::Size(w, h), CV_32F);
   for (unsigned int i = 0; i < h; i++)
   {
     for (unsigned int j = 0; j < w; j++)
     {
       output.at<float>(i, j) = iter.Get(); // TBD: remove ".at" handler, instead of i & j loops, just use 'iter' and 'iter.GetIndex()' for looping
       ++iter;
     }

   }
   //  output = out;
 }

 ImageTypeFloat3D::Pointer convertOpencvMattoITK(cv::Mat input)
 {
   auto new_image = ImageTypeFloat3D::New();
   ImageTypeFloat3D::RegionType region;
   ImageTypeFloat3D::RegionType::SizeType size;
   ImageTypeFloat3D::RegionType::IndexType start;
   ImageTypeFloat3D::PointType origin;
   origin.Fill(0);
   ImageTypeFloat3D::SpacingType spacing;
   size.Fill(1);
   size[0] = input.cols;
   size[1] = input.rows;
   start.Fill(0);
   spacing.Fill(1);
   region.SetSize(size);
   region.SetIndex(start);
   new_image->SetRegions(region);
   new_image->SetRequestedRegion(region);
   new_image->SetBufferedRegion(region);
   new_image->Allocate();
   new_image->FillBuffer(0);
   new_image->SetOrigin(origin);
   new_image->SetSpacing(spacing);
   IteratorType iter(new_image, new_image->GetLargestPossibleRegion());
   iter.GoToBegin();
   for (int i = 0; i < input.rows; i++)
   {
     for (int j = 0; j < input.cols; j++)
     {
       iter.Set(input.at<float>(i, j));
       ++iter;
     }

   }
   return new_image;
 }


};

