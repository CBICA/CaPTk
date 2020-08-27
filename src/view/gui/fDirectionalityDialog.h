///////////////////////////////////////////////////////////////////////////////////////
// fDirectionalityDialog.h
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

#ifndef _fDirectionalityDialog_h_
#define _fDirectionalityDialog_h_


//#include "CAPTk.h"
#include "ui_fDirectionalityDialog.h"

/**
\class fDirectionalityDialog

\brief This class controls the elements in the recurrence dialog
*/
class fDirectionalityDialog : public QDialog, private Ui::fDirectionalityDialog
{
  Q_OBJECT

public:
  fDirectionalityDialog();
  ~fDirectionalityDialog();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void OpenROI1Image();
  void OpenROI2Image();
  void SelectOutputFolder();

signals:
  void RunDirectionalityEstimator(const std::string roi1Image, const std::string roi2Image, const std::string outputDir);
};


#endif
