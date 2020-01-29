///////////////////////////////////////////////////////////////////////////////////////
// fRecurrenceDialog.h
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

#ifndef _fRecurrenceDialog_h_
#define _fRecurrenceDialog_h_


//#include "CAPTk.h"
#include "ui_fRecurrenceDialog.h"

#define SUBJECT_CLASSIFICATION 0
#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fRecurrenceDialog

\brief This class controls the elements in the recurrence dialog
*/
class fRecurrenceDialog : public QDialog, private Ui::fRecurrenceDialog
{
  Q_OBJECT

public:
  fRecurrenceDialog();
  ~fRecurrenceDialog();
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
  void OpenSVMModelFile1(); 
  void OpenSVMModelFile2();
  void OpenExistingMasksDirectory();
  void SelectOutputDirectory();
  void ExistingClassificationRadioButtonChecked();
  void LoadedClassificationRadioButtonChecked();
  void NewModelRadioButtonChecked();
  void OpenTestSubjectsDirectory();
  void CurrentSubjectRadioButtonChecked();
  void CheckForDisclaimer();

signals:
  void SubjectBasedRecurrenceEstimate(std::string outputdirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistanceData);
  void SubjectBasedExistingRecurrenceEstimate(std::string outputdirectory, std::string modeldirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistanceData);
  void ExistingModelBasedRecurrenceEstimate(const std::string &modeldirectory, const std::string &inputdirectory, const std::string &outputdirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistanceData);
  void TrainNewModel(std::string directory, std::string outputdirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistanceData);
};


#endif





