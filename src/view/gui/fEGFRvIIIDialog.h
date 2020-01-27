///////////////////////////////////////////////////////////////////////////////////////
// fEGFRvIIIPredictor.h
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
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fEGFRvIIIPredictor_h_
#define _fEGFRvIIIPredictor_h_


//#include "CAPTk.h"
#include "ui_fEGFRvIIIDialog.h"

#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fEGFRvIIIPredictor

\brief This class controls the elements in the recurrence dialog
*/
class fEGFRvIIIPredictor : public QDialog, private Ui::fEGFRvIIIPredictor
{
  Q_OBJECT

public:
  fEGFRvIIIPredictor();
  ~fEGFRvIIIPredictor();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }
  void SetCurrentLoggerPath(std::string loggerPath)
  {
    mLoggerFile = loggerPath;
  }

  void SetTrainedModelLink(std::string link)
  {
    m_trainedModelLink = link;
  }
  std::string m_trainedModelLink;

  QString mInputPathName;
  std::string mLoggerFile;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenSVMModelFile();
  void OpenExistingMasksDirectory();
  void SelectOutputDirectory();
  void ExistingClassificationRadioButtonChecked();
  void NewModelRadioButtonChecked();
  void OpenTestSubjectsDirectory();
  void CheckForDisclaimer();

signals:
  void PrepareNewEGFRvIIIPredictionModel(const std::string inputdirectory, const std::string outputdirectory);
  void EGFRvIIIPredictionOnExistingModel(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);
  void SubjectBasedSurvivalEstimate(const std::string output, const std::string model, double age);
};


#endif
