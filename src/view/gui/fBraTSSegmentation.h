///////////////////////////////////////////////////////////////////////////////////////
// fBraTSSegmentation.h
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
#pragma once


//#include "CAPTk.h"
#include "ui_fBraTSSegmentation.h"

#include <map>

/**
\class fBraTSSegmentation

\brief This class controls the elements in the recurrence dialog
*/
class fBraTSSegmentation : public QDialog, private Ui::fBraTSSegmentation
{
  Q_OBJECT

public:
  fBraTSSegmentation();
  ~fBraTSSegmentation();

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    //outputImageName->setText(mInputPathName + "/output");
  }

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void OpenInputT1CEImage();
  void OpenInputT1Image();
  void OpenInputT2Image();
  void OpenInputFLImage();
  void SelectOutputDirectory();

signals:
  void RunBraTSPipeline(
    const std::string inputT1CEImageFile, 
    const std::string inputT1ImageFile,
    const std::string inputT2ImageFile,
    const std::string inputFLImageFile,
    const std::string outputImageFile
  );
};

