///////////////////////////////////////////////////////////////////////////////////////
// fSurvivalPredictor.h
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

#ifndef _fSurvivalPredictor_h_
#define _fSurvivalPredictor_h_


#include "CAPTk.h"
#include "ui_fSurvivalDialog.h"

#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fSurvivalPredictor

\brief This class controls the elements in the recurrence dialog
*/
class fSurvivalPredictor : public QDialog, private Ui::fSurvivalPredictor
{
  Q_OBJECT

public:
  fSurvivalPredictor();
  ~fSurvivalPredictor();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenSVMModelFile();
  void OpenExistingMasksDirectory();
  void SelectOutputDirectory();
  void ExistingClassificationRadioButtonChecked();
  void NewModelRadioButtonChecked();
  void OpenTestSubjectsDirectory();

signals:
  void PrepareNewSurvivalPredictionModel(const std::string inputdirectory, const std::string outputdirectory);
  void SurvivalPredictionOnExistingModel(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);
  void SubjectBasedSurvivalEstimate(const std::string output, const std::string model, double age);
};


#endif
