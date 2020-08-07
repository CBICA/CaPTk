///////////////////////////////////////////////////////////////////////////////////////
// fFeaturePanel.cpp
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
// Contact details: software@cbica.upenn.edu
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////
#include "fFeaturePanel.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "fMainWindow.h"
#include "FeatureExtraction/src/FeatureExtraction.h"
#include "exception"

#include "itkExtractImageFilter.h"


fFeaturePanel::fFeaturePanel(QWidget * parent) : QWidget(parent)
{
  m_listener = NULL;
  m_featureDialog = NULL;
  setupUi(this);

  m_featureFiles = cbica::stringSplit(std::string(FeatureDefinitions), " ");
  m_featureFiles.erase(m_featureFiles.begin()); // because an empty string always gets appended first

  connect(m_btnCompute, SIGNAL(clicked()), this, SLOT(onComputeButtonClicked()));
  connect(m_btnAdvanced, SIGNAL(clicked()), this, SLOT(advancedButtonClicked()));
  connect(m_btnBrowseSaveFile, SIGNAL(clicked()), this, SLOT(browseOutputFileName()));
  connect(m_cmbFeatureType, SIGNAL(currentIndexChanged(int)), this, SLOT(featureTypeChanged(int)));
  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  //csv_format->setChecked(true);
  m_verticalConcat->setChecked(true);
  radio1->setChecked(true);
  loadFeatureFiles();
  featureTypeChanged(0);

}

void fFeaturePanel::helpClicked()
{
  emit helpClicked_FeaUsage("ht_FeatureExtraction.html");
}
std::map< std::string, bool > fFeaturePanel::getEnabledFeatures()
{
  std::map< std::string, bool > result;
  int type = m_cmbFeatureType->currentIndex();
  if (type >= m_FeatureMaps.size())
  {
    return std::map< std::string, bool >();
  }

  for (auto chk : m_featureCheckBoxMap)
  {
    if (chk.first == "Intensity")
    {
      chk.second->setChecked(true);
      chk.second->setEnabled(false);

      result[chk.first] = chk.second->isChecked();
    }
    //else if (chk.first == "LBP2D")
    //{
    //  chk.second->setChecked(false);
    //  chk.second->setEnabled(false);
    //}

    if (!m_FeatureMaps[type][chk.first].empty())
    {
      result[chk.first] = chk.second->isChecked();
    }
  }
  return result;
}
void fFeaturePanel::featureTypeChanged(int type)
{
  if (type >= m_FeatureMaps.size())
  {
    return;
  }
  if (type >= m_featureFiles.size())
  {
    return;
  }

  auto currentMap = m_FeatureMaps[type];
  if (!currentMap.empty())
  {
    for (auto chk : m_featureCheckBoxMap)
    {
      std::string name = chk.first;
      //ckh.second->setChecked(!currentMap[ckh.first].empty());
      //ckh.second->setEnabled(!currentMap[ckh.first].empty());
      if (currentMap[name].empty())
      {
        chk.second->setChecked(false);
        chk.second->setEnabled(false);
      }
      else
      {
        if (name == "Intensity")
        {
          chk.second->setChecked(true);
          chk.second->setEnabled(false);
        }
        else
        {
          chk.second->setChecked(true);
          chk.second->setEnabled(true);
        }
      }
    }
  }

}
void fFeaturePanel::browseOutputFileName()
{
  QFileDialog dialog(this, tr("Save Features"), loggerFolder.c_str(), tr("Feature File (*.csv *.xml)"));
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.selectFile("features.csv");

  int ret = dialog.exec();
  if (ret == QDialog::Accepted)
  {
    QString fileName = dialog.selectedFiles()[0];
    m_txtSaveFileName->setText(fileName);
  }
}
void fFeaturePanel::ComputeFunctionality()
{

  int  type = m_cmbFeatureType->currentIndex();
  m_cmbFeatureType->setEnabled(false);
  computeFeature(type);
  m_cmbFeatureType->setEnabled(true);
}
void fFeaturePanel::CancelFunctionality()
{
  //  features_extraction.cancel();

}
void fFeaturePanel::onComputeButtonClicked()
{
  ComputeFunctionality();
}
/*When compute button clicked this function is called*/
void fFeaturePanel::computeFeature(int type)
{

  std::vector < ImageTypeFloat3D::Pointer> images;
  ImageTypeFloat3D::Pointer maskImg;
  std::vector< std::string >imagepaths;
  std::vector<std::string> modality;


  if (m_listener != NULL)
  {
    if (radio1->isChecked())
    {
      images.push_back(
        ((fMainWindow*)m_listener)->mSlicerManagers[
          ((fMainWindow*)m_listener)->mCurrentPickedImageIndex]->GetITKImage()
            );
      modality.push_back(
        CAPTK::ImageModalityString[
          ((fMainWindow*)m_listener)->mSlicerManagers[
            ((fMainWindow*)m_listener)->mCurrentPickedImageIndex]->mImageSubType]
      );
    }
    else
    {
      images = ((fMainWindow*)m_listener)->getLodedImages(imagepaths, modality);
      if (imagepaths.size() < 1)
      {
        ShowErrorMessage("No valid images selected!", this);
        return;
      }
      if (images.size() != imagepaths.size())
      {
        ShowErrorMessage("No valid images selected!", this);
        return;
      }
    }
    //((fMainWindow*)m_listener)->mCurrentPickedImageIndex;
    if (((fMainWindow*)m_listener)->isMaskDefined())
    {
      maskImg = ((fMainWindow*)m_listener)->getMaskImage();
    }
    else
    {
      ShowErrorMessage("No valid ROI defined!", this);
      return;
    }
  }
  if (images.size() == 0)
  {
    ShowErrorMessage("No valid images selected!", this);
    return;
  }
  else
  {
    ((fMainWindow*)m_listener)->updateProgress(5, "Initializing Feature Extraction Module");
  }

  std::string featureFileName = m_txtSaveFileName->text().toStdString();

  /*  if (cbica::getFilenameExtension(featureFileName).empty())
  {
  if (cbica::fileExists(featureFileName + ".csv"))
  {
  ShowErrorMessage("File already exists. Please give a new filename");
  return;
  }
  if (cbica::fileExists(featureFileName + ".xml"))
  {
  ShowErrorMessage("File already exists. Please give a new filename");
  return;
  }
  }
  else*/ if (cbica::fileExists(featureFileName))
  {
    ShowErrorMessage("File already exists. Please give a new filename");
    ((fMainWindow*)m_listener)->updateProgress(0, "");
    return;
  }

  std::string roi_given = m_roi->text().toStdString();
  std::string roi_label_given = m_roi_label->text().toStdString();

  auto featureMap = m_FeatureMaps[type];
  auto selectedFeatures = getEnabledFeatures();

  // check for resampling rate
  for (auto &currentFeature : featureMap) // iterating over the feature families, i.e., GLCM, Intensity, ...
  {
    if (currentFeature.first.find("Generic") != std::string::npos)
    {
      for (size_t j = 0; j < currentFeature.second.size(); j++)
      {
        for (auto &currentFeature_Parameter : currentFeature.second[j]) // each parameter within the feature family
        {
          if (currentFeature_Parameter.first.find("Value"))
          {
            if (std::atof(currentFeature_Parameter.second.c_str()) == 3.0)
            {
              ShowMessage(
                "The resampling rate is set to 3.0 for computation efficiency, please reduce it to increase accuracy of extracted features\n\nChange it under 'Customization -> Advanced ... -> Generic -> Resampling Rate'",
                this, "Warning");
            }
          }
        }
      }
    }
  }

  ((fMainWindow*)m_listener)->updateProgress(30, "Calculating and exporting features", images.size());

  if (images[0]->GetLargestPossibleRegion().GetSize()[2] == 1)
  {
    // this is actually a 2D image which has been loaded as a 3D image with a single slize in z-direction

    cbica::Logging(loggerFile, "2D Image detected, doing conversion and then passing into FE module");
    using ImageTypeFloat2D = itk::Image< float, 2 >;
    std::vector< ImageTypeFloat2D::Pointer > images_2d(images.size());

    ImageTypeFloat3D::IndexType regionIndex;
    regionIndex.Fill(0);
    auto regionSize = images[0]->GetLargestPossibleRegion().GetSize();
    regionSize[2] = 0; // only 2D image is needed
    for (size_t i = 0; i < images.size(); i++)
    {
      auto extractor = itk::ExtractImageFilter< ImageTypeFloat3D, ImageTypeFloat2D >::New();
      ImageTypeFloat3D::RegionType desiredRegion(regionIndex, regionSize);
      extractor->SetExtractionRegion(desiredRegion);
      extractor->SetInput(images[i]);
      extractor->SetDirectionCollapseToIdentity();
      extractor->Update();
      images_2d[i] = extractor->GetOutput();
      images_2d[i]->DisconnectPipeline();
    }

    auto extractor = itk::ExtractImageFilter< ImageTypeFloat3D, ImageTypeFloat2D >::New();
    ImageTypeFloat3D::RegionType desiredRegion(regionIndex, regionSize);
    extractor->SetExtractionRegion(desiredRegion);
    extractor->SetInput(maskImg);
    extractor->SetDirectionCollapseToIdentity();
    extractor->Update();
    auto mask2D = extractor->GetOutput();
    //mask2D->DisconnectPipeline();

    FeatureExtraction<  ImageTypeFloat2D >  featureextractor;
    featureextractor.SetPatientID(m_patientID_label->text().toStdString());
    featureextractor.SetInputImages(images_2d, modality);
    featureextractor.SetMaskImage(mask2D);
    featureextractor.SetVerticallyConcatenatedOutput(true);
    featureextractor.SetRequestedFeatures(featureMap, selectedFeatures);
    featureextractor.SetSelectedROIsAndLabels(roi_given, roi_label_given);
    featureextractor.SetOutputFilename(m_txtSaveFileName->text().toStdString());
    try
    {
      featureextractor.Update();
    }
    catch (const std::exception& e)
    {
      cbica::Logging(loggerFile, e.what());
      ShowErrorMessage(e.what());
    }  
  }
  else
  {
    FeatureExtraction<  ImageTypeFloat3D >  featureextractor;
    featureextractor.SetPatientID(m_patientID_label->text().toStdString());
    featureextractor.SetInputImages(images, modality);
    featureextractor.SetMaskImage(maskImg);
    featureextractor.SetVerticallyConcatenatedOutput(m_verticalConcat->isChecked());
    featureextractor.SetRequestedFeatures(featureMap, selectedFeatures);
    featureextractor.SetSelectedROIsAndLabels(roi_given, roi_label_given);
    featureextractor.SetOutputFilename(m_txtSaveFileName->text().toStdString());
    try
    {
      featureextractor.Update();
    }
    catch (const std::exception& e)
    {
      cbica::Logging(loggerFile, e.what());
      ShowErrorMessage(e.what());
    }
  }

  if (/*!featurevec.empty() && */!featureFileName.empty())
  {

    //if (csv_format->checkState() == Qt::Checked)
    //{
    //  std::string filetype = "csv";
    //  //      features_extraction.writeFeatureList(featureFileName, featurevec, filetype);
    //}

    //if ((csv_format->checkState() != Qt::Checked) && (xml_format->checkState() != Qt::Checked))
    //{
    //  ShowErrorMessage("No format type selected; writing as XML");
    //  csv_format->isEnabled();
    //}

    //if (xml_format->checkState() == Qt::Checked)
    //{
    //  std::string filetype = "xml";
    //  //      features_extraction.writeFeatureList(featureFileName, featurevec, filetype);
    //}

    ((fMainWindow*)m_listener)->updateProgress(100, "Feature export complete!");
  }
  else
  {
    //if (featurevec.empty())
    //{
    //  ShowErrorMessage("File to write is empty, please check");
    //}
  }

  m_btnCompute->setEnabled(true);
  return;
}

void fFeaturePanel::advancedButtonClicked()
{
  // This is not finalized code - need to save the feature tuple and need to be used with compute function 
  int selectedFeatureType = m_cmbFeatureType->currentIndex();
  auto featureMap = m_FeatureMaps[selectedFeatureType];
  if (m_featureDialog)
  {
    delete m_featureDialog;
  }
  m_featureDialog = new FeatureDialogTree(this);
  m_featureDialog->setModal(true);
  if (featureMap.empty())
  {
    m_featureDialog->showMessage("Empty feature map: ");
  }
  else
  {
    m_featureDialog->addFeatureMap(featureMap, getEnabledFeatures());
  }
  m_featureDialog->exec();
  m_FeatureMaps[selectedFeatureType] = m_featureDialog->getFeatureMap();

}