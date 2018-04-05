///////////////////////////////////////////////////////////////////////////////////////
// fFeaturePanel.h
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
// License Agreement: http://www.med.upenn.edu/sbia/software/license.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fFeaturePanel_h_
#define _fFeaturePanel_h_


#include "CAPTk.h"
#include "Landmarks.h"
#include "ui_fFeaturePanel.h"
#include "fFeatureDialog.h"
#include "FeatureExtraction.h"
/**
\class fFeaturePanel

\brief This class controls the elements in the Feature panel of the tab
*/
class fFeaturePanel : public QWidget, private Ui::fFeaturePanel
{
  Q_OBJECT

public:
  //! Constructor
  fFeaturePanel(QWidget * parent = 0);

  //! Destructor
  ~fFeaturePanel() {}
  // void  writeFeatureList(std::string Filename, std::vector< std::tuple<std::string, std::string, float>>featurevec);

  template <class TImageType = itk::Image< float, 3 >>
  typename TImageType::Pointer get_selected_mask(typename TImageType::Pointer total_mask, int roi);

  void setListner(void* lst)
  {
    m_listener = lst;

  }

  //! Sets the same temporary folder everywhere
  void setTempFolderLocation(const std::string& input_tempFolder)
  {
    tempFolderLocation = input_tempFolder;
  }


signals:
  void m_btnComputeClicked();
  void helpClicked_FeaUsage(std::string);

  public slots :
  void browseOutputFileName();
  void ComputeFunctionality();
  void CancelFunctionality();
  void onComputeButtonClicked();
  void computeFeature(int type);
  void featureTypeChanged(int type);
  void advancedButtonClicked();
  void helpClicked();
  map<string, bool> getEnabledFeatures();

private:
  //  FeatureExtraction features_extraction; // TBD change i to m_ 
  FeatureDialogTree* m_featureDialog;
  vector < map<string, vector<map<string, string> > > > m_FeatureMaps;// TBD typedef the map value so the code doenst look ugly 
  void* m_listener;//TBD this is a bad design (because of time pressure): needs to be changed RK
private:
  void loadFeatureFiles() // TBD move to cpp 
  {
    //auto names = FeatureExtraction::getFeatureMapFiles();

    std::string dataFeatureDir = cbica::normPath(QApplication::applicationDirPath().toStdString() + "/../data/features/");

    //std::string defaultFeatureFile = dataFeatureDir + "/1_params_default.csv";
    if (!cbica::isFile(dataFeatureDir + "/1_params_default.csv"))
    {
      dataFeatureDir = cbica::normPath(QApplication::applicationDirPath().toStdString() + "/../../data/features/");
      //defaultFeatureFile = dataFeatureDir + "/1_params_default.csv";
    }

    auto filesInDir = cbica::filesInDirectory(dataFeatureDir);
    for (size_t i = 0; i < filesInDir.size(); i++)
    {
      auto featureMap = FeatureParser(filesInDir[i]).getFeatureMap();
      m_FeatureMaps.push_back(featureMap);
    }
  }

};


#endif
