///////////////////////////////////////////////////////////////////////////////////////
// fFeaturePanel.cpp
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
// Contact details: software@cbica.upenn.edu
// License Agreement: http://www.med.upenn.edu/sbia/software/license.html
///////////////////////////////////////////////////////////////////////////////////////
#include "fFeaturePanel.h"
#include "CAPTk.h"
#include "fMainWindow.h"
#include "FeatureExtraction.h"
#include "exception"



fFeaturePanel::fFeaturePanel(QWidget * parent) : QWidget(parent)
{
  m_listener = NULL;
  m_featureDialog = NULL;
  setupUi(this);
  connect(m_btnCompute, SIGNAL(clicked()), this, SLOT(onComputeButtonClicked()));
  connect(m_btnAdvanced, SIGNAL(clicked()), this, SLOT(advancedButtonClicked()));
  connect(m_btnBrowseSaveFile, SIGNAL(clicked()), this, SLOT(browseOutputFileName()));
  connect(m_cmbFeatureType, SIGNAL(currentIndexChanged(int)), this, SLOT(featureTypeChanged(int)));
  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  csv_format->setChecked(true);
  radio1->setChecked(true);
  loadFeatureFiles();
  featureTypeChanged(0);

  //m_NeighbGrayTone->hide();
  //m_Runlenghtfeature->hide();
  //m_GraylevelSizeZone->hide();

  //m_NeighbGrayTone->setChecked(false);
  //m_Runlenghtfeature->setChecked(false);
  //m_GraylevelSizeZone->setChecked(false);
}

void fFeaturePanel::helpClicked()
{
  emit helpClicked_FeaUsage("ht_FeatureExtraction.html");
}
map<string, bool> fFeaturePanel::getEnabledFeatures()
{
  map<string, bool> result;
  int type = m_cmbFeatureType->currentIndex();
  if (type >= m_FeatureMaps.size())
  {
    return map<string, bool>();
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
  auto currentMap = m_FeatureMaps[type];
  if (!currentMap.empty())
  {
    for (auto chk : m_featureCheckBoxMap)
    {
      string name = chk.first;
      //ckh.second->setChecked(!currentMap[ckh.first].empty());
      //ckh.second->setEnabled(!currentMap[ckh.first].empty());
      if (currentMap[name].empty())
      {
        chk.second->setChecked(false);
        chk.second->setEnabled(false);
      }
      else
      {
        if (chk.first == "Intensity")
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
    images = ((fMainWindow*)m_listener)->getLodedImages(imagepaths, modality);
    if (imagepaths.size() <1)
    {
      ShowErrorMessage("No valid images selected!", this);
      return;
    }
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
  if (images.size() != imagepaths.size() || images.size() == 0)
  {
    ShowErrorMessage("No valid images selected!", this);
    return;
  }
  else
  {
    ((fMainWindow*)m_listener)->updateProgress(5, "Computing features!");
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
    return;
  }

  std::string roi_given = m_roi->text().toStdString();
  std::string roi_label_given = m_roi_label->text().toStdString();

  auto featureMap = m_FeatureMaps[type];
  auto selectedFeatures = getEnabledFeatures();

  std::string feature_selected;




  typedef itk::VectorContainer< unsigned char, OffsetType > OffsetVector;
  typedef OffsetVector::Pointer  OffsetVectorPointer;
  //OffsetType offset;
  FeatureExtraction< ImageTypeFloat3D, float>  featureextractor;
  featureextractor.SetInputImages(images, modality);
  featureextractor.SetMaskImage(maskImg);

  std::string featureFile = cbica::normPath(QApplication::applicationDirPath().toStdString() + "/../data/features/params_default.csv");
  //std::string dataDir = "E:/Projects/newcaptk/data/features/params_default.csv";
  if (!cbica::isFile(featureFile))
  {
    featureFile = cbica::normPath(QApplication::applicationDirPath().toStdString() + "/../../data/features/params_default.csv");
    //dataDir = "E:/Projects/newcaptk/data/features/params_default.csv";
  }

  featureextractor.SetRequestedFetures(featureMap, selectedFeatures);
  //featureextractor.SetRequestedFetures(featureFile, selectedFeatures);
  featureextractor.SetSelectedLabels(roi_given, roi_label_given);
  featureextractor.Setouputfilename(m_txtSaveFileName->text().toStdString());
  featureextractor.Update();

  if (/*!featurevec.empty() && */!featureFileName.empty())
  {

    if (csv_format->checkState() == Qt::Checked)
    {
      std::string filetype = "csv";
      //      features_extraction.writeFeatureList(featureFileName, featurevec, filetype);
    }

    if ((csv_format->checkState() != Qt::Checked) && (xml_format->checkState() != Qt::Checked))
    {
      ShowErrorMessage("No format type selected; writing as XML");
      csv_format->isEnabled();
    }

    if (xml_format->checkState() == Qt::Checked)
    {
      std::string filetype = "xml";
      //      features_extraction.writeFeatureList(featureFileName, featurevec, filetype);
    }


    ((fMainWindow*)m_listener)->updateProgress(70, "Calculating and exporting features", images.size());

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