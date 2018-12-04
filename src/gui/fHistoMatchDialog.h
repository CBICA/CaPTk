///////////////////////////////////////////////////////////////////////////////////////
// fHistoMatcher.h
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

#ifndef _fHistoMatchDialog_h_
#define _fHistoMatchDialog_h_


//#include "CAPTk.h"
#include "ui_fHistoMatchDialog.h"

/**
\class fHistoMatcher

\brief This class controls the elements in the recurrence dialog
*/
class fHistoMatcher : public QDialog, private Ui::fHistoMatcher
{
  Q_OBJECT

public:
  fHistoMatcher();
  ~fHistoMatcher();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName + "histoMatch_output.nii.gz");
  }

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void OpenReferenceImage();
  void OpenInputImage();
  void SelectOutputImage();

signals:
  void RunHistogramMatching(const std::string referenceFile, const std::string inputImageFile, const std::string outputImageFile);
};


#endif
