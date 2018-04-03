///////////////////////////////////////////////////////////////////////////////////////
// fRegistrationDialog.h
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

#ifndef _fRegistrationDialog_h_
#define _fRegistrationDialog_h_


#include "CAPTk.h"
#include "ui_fRegistrationDialog.h"

#define SUBJECT_CLASSIFICATION 0
#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fRegistrationDialog

\brief This class controls the elements in the registration dialog
*/
class fRegistrationDialog : public QDialog, private Ui::fRegistrationDialog
{
  Q_OBJECT

public:
  fRegistrationDialog();
  ~fRegistrationDialog();

  QString mInputPathName;

  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void SelectFixedFile();
  void SelectMovingFile1();
  void SelectMovingFile2();
  void SelectMovingFile3();
  void SelectMovingFile4();
  void SelectMovingFile5();
  void SelectMovingOutputFile1();
  void SelectMovingOutputFile2();
  void SelectMovingOutputFile3();
  void SelectMovingOutputFile4();
  void SelectMovingOutputFile5();

signals:
  void Registrationsignal(std::string fixedfilename, std::vector<std::string> inputfilenames, std::vector<std::string> outputfilenames);
};


#endif





