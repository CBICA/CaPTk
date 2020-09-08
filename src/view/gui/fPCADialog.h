///////////////////////////////////////////////////////////////////////////////////////
// fPCADialog.h
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

#ifndef _fPCADialog_h_
#define _fPCADialog_h_


//#include "CAPTk.h"
#include "ui_fPCADialog.h"

#define SUBJECT_CLASSIFICATION 0
#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fPCADialog

\brief This class controls the elements in the recurrence dialog
*/
class fPCADialog : public QDialog, private Ui::fPCADialog
{
  Q_OBJECT

public:
  fPCADialog();
  ~fPCADialog();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }
  void SetCurrentLoggerPath(std::string loggerPath)
  {
    mLoggerFile = loggerPath;
  }

  QString mInputPathName;
  std::string mLoggerFile;

public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenSVMModelFile1();
  void OpenSVMModelFile2();
  void OpenExistingMasksDirectory();
  void SelectOutputDirectory();
  void ExistingClassificationRadioButtonChecked();
  void LoadedClassificationRadioButtonChecked();
  void NewModelRadioButtonChecked();
  void OpenTestSubjectsDirectory();
  void CurrentSubjectRadioButtonChecked();
  //void CheckForDisclaimer();

signals:
  void ExistingModelBasedPCAEstimate(const std::string &modeldirectory, const std::string &inputdirectory, const std::string &outputdirectory);
  void TrainNewPCAModel(std::string directory, std::string outputdirectory);
};


#endif





