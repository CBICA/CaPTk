///////////////////////////////////////////////////////////////////////////////////////
// fTexturePipelineDialog.h
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

#ifndef _fTexturePipelineDialog_h_
#define _fTexturePipelineDialog_h_


//#include "CAPTk.h"
#include "ui_fTexturePipelineDialog.h"

/**
\class fTexturePipelineDialog

\brief This class controls the elements in the DICOM converter
*/
class fTexturePipelineDialog : public QDialog, private Ui::fTexturePipelineDialog
{
  Q_OBJECT

public:

  fTexturePipelineDialog();
  ~fTexturePipelineDialog();

  QString m_exe, m_dataDir, m_modelDir; // contains full path and exe name of dcm2nii

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }
  
  QString mInputPathName;
  
public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void SelectOutputDirectory();

signals:
  void RunTextureFeaturePipeline(const std::string outputDirectory);
};


#endif
