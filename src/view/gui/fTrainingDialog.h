///////////////////////////////////////////////////////////////////////////////////////
// fTrainingSimulator.h
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

#ifndef _fTrainingDialog_h_
#define _fTrainingDialog_h_


#include "CaPTkEnums.h"
#include "ui_fTrainingDialog.h"
#include "TrainingModuleParameters.h"

/**
\class fTrainingSimulator

\brief This class controls the elements in the recurrence dialog
*/
class fTrainingSimulator : public QDialog, private Ui::fTrainingSimulator
{
  Q_OBJECT

public:
  fTrainingSimulator();
  ~fTrainingSimulator();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputFeaturesName = inputPath;
    outputDirectoryName->setText(mInputFeaturesName);
  }

  QString mModelDirectoryName;
  QString mInputFeaturesName;
  QString mInputTargetName;
  QString mInputBValName;
  QString mInputBVecName;
  QString mOutputPathName;



public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenInputImage();
  void OpenInputMaskImage();
  void SelectOutputImage();
  void SelectSplitModelDirectory();
  void CrossValidationRadioButtonChecked();
  void SplitTrainRadioButtonChecked();
  void SplitTestRadioButtonChecked();
  void OptimizationToggled(bool);

signals:
  void RunTrainingSimulation(const TrainingModuleParameters params);
};


#endif
