///////////////////////////////////////////////////////////////////////////////////////
// fMolecularSubtypePredictor.h
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

#ifndef _fMolecularSubtypePredictor_h_
#define _fMolecularSubtypePredictor_h_


//#include "CAPTk.h"
#include "ui_fMolecularSubtypeDialog.h"

#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fMolecularSubtypePredictor

\brief This class controls the elements in the recurrence dialog
*/
class fMolecularSubtypePredictor : public QDialog, private Ui::fMolecularSubtypePredictor
{
  Q_OBJECT

public:
  fMolecularSubtypePredictor();
  ~fMolecularSubtypePredictor();
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
  void PrepareNewMolecularSubtypePredictionModel(const std::string inputdirectory, const std::string outputdirectory);
  void MolecularSubtypePredictionOnExistingModel(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);
  void SubjectBasedSurvivalEstimate(const std::string output, const std::string model, double age);
};


#endif
