///////////////////////////////////////////////////////////////////////////////////////
// fSegmentationPanel.h
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/cbica/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fSegmentationPanel_h_
#define _fSegmentationPanel_h_


#include "CAPTk.h"
#include "Landmarks.h"
#include "ui_fSegmentationPanel.h"
#include "fFeatureDialog.h"
#include "FeatureExtraction.h"
/**
\class fSegmentationPanel

\brief This class controls the elements in the Feature panel of the tab
*/
class fSegmentationPanel : public QWidget, private Ui::fSegmentationPanel
{
  Q_OBJECT

public:
  //! Constructor
  fSegmentationPanel(QWidget * parent = 0) {};

  //! Destructor
  ~fSegmentationPanel() {};
  // void  writeFeatureList(std::string Filename, std::vector< std::tuple<std::string, std::string, float>>featurevec);

  void setListner(void* lst)
  {
    //m_listener = lst;
  }

  //! Sets the same temporary folder everywhere
  void setTempFolderLocation(const std::string& input_tempFolder)
  {
    m_tempFolderLocation = input_tempFolder;
  }


signals:
  //void m_btnComputeClicked();
  //void helpClicked_FeaUsage(std::string);

public slots:
 /* void browseOutputFileName();
  void ComputeFunctionality();
  void CancelFunctionality();
  void onComputeButtonClicked();
  void computeFeature(int type);
  void featureTypeChanged(int type);
  void advancedButtonClicked();
  void helpClicked();
  std::map< std::string, bool > getEnabledFeatures();*/

private:
  ////  FeatureExtraction features_extraction; // TBD change i to m_ 
  //FeatureDialogTree* m_featureDialog;
  //std::vector< std::map< std::string, std::vector< std::map< std::string, std::string> > > > m_FeatureMaps;// TBD typedef the map value so the code doenst look ugly 
  //void* m_listener;//TBD this is a bad design (because of time pressure): needs to be changed RK
  //void loadFeatureFiles() // TBD move to cpp 
  //{
  //  //auto names = FeatureExtraction::getFeatureMapFiles();

  //  std::string dataFeatureDir = cbica::normPath(getCaPTkDataDir() + "/features");

  //  //std::string defaultFeatureFile = dataFeatureDir + "/1_params_default.csv";
  //  if (!cbica::isFile(dataFeatureDir + "/" + m_featureFiles[0]))
  //  {
  //    dataFeatureDir = cbica::normPath(captk_currentApplicationPath + "/../../data/features/");
  //    //defaultFeatureFile = dataFeatureDir + "/1_params_default.csv";
  //  }

  //  auto filesInDir = cbica::filesInDirectory(dataFeatureDir);
  //  if (filesInDir.size() != m_featureFiles.size())
  //  {
  //    cbica::Logging(loggerFile, "Feature file number mismatch");
  //    return;
  //  }
  //  for (size_t i = 0; i < filesInDir.size(); i++)
  //  {
  //    auto featureMap = FeatureParser(filesInDir[i]).getFeatureMap();
  //    m_FeatureMaps.push_back(featureMap);
  //  }
  //}

};


#endif
