///////////////////////////////////////////////////////////////////////////////////////
// fPreprocessingDialog.h
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

#ifndef _fPreprocessingDialog_h_
#define _fPreprocessingDialog_h_


//#include "CAPTk.h"
#include "ui_fPreprocessingDialog.h"

#define SUBJECT_CLASSIFICATION 0
#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fPreprocessingDialog

\brief This class controls the elements in the pre-processing dialog
*/
class fPreprocessingDialog : public QDialog, private Ui::fPreprocessingDialog
{
  Q_OBJECT

public:
  fPreprocessingDialog();
  ~fPreprocessingDialog();
  int mode;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenExistingMasksDirectory();
  void OpenSVMModelFile();
  void OpenTestSubjectsDirectory();
  void SelectOutputDirectory();
  void ExistingClassificationRadioButtonChecked();
  void NewModelRadioButtonChecked();


signals:

};


#endif





