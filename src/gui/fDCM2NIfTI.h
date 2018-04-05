///////////////////////////////////////////////////////////////////////////////////////
// fDCM2NIfTI.h
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

#ifndef _fDCM2NIfTI_h_
#define _fDCM2NIfTI_h_


#include "CAPTk.h"
#include "ui_fDCM2NIfTI.h"

/**
\class fDCM2NIfTI

\brief This class controls the elements in the DICOM converter
*/
class fDCM2NIfTIConverter : public QDialog, private Ui::fDCM2NIfTIConverter
{
  Q_OBJECT

public:
  fDCM2NIfTIConverter();
  ~fDCM2NIfTIConverter();

  QString m_exe; // contains full path and exe name of dcm2nii

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName + "dcm2nii.nii.gz");
  }

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void OpenInputImage();
  void SelectOutputImage();

signals:
  void RunDICOMConverter(const std::string firstImageInSeries, const std::string outputImage);
};


#endif
