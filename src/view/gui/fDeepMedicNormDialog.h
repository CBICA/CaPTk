///////////////////////////////////////////////////////////////////////////////////////
// fDeepMedicNormalizer.h
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

#ifndef _fDeepMedicNormDialog_h_
#define _fDeepMedicNormDialog_h_


//#include "CAPTk.h"
#include "ui_fDeepMedicNormDialog.h"

/**
\class fDeepMedicNormalizer

\brief This class controls the elements in the recurrence dialog
*/
class fDeepMedicNormalizer : public QDialog, private Ui::fDeepMedicNormalizer
{
  Q_OBJECT

public:
  fDeepMedicNormalizer();
  ~fDeepMedicNormalizer();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName + "deepMedicNorm_output.nii.gz");
  }

  QString mInputPathName;

  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void OpenMaskImage();
  void OpenInputImage();
  void SelectOutputImage();

signals:
  void RunDeepMedicNormalizer(const std::string inputImageFile, const std::string maskFile, const std::string outputImageFile,
    const std::string quantLower, const std::string quantUpper,
    const std::string cutoffLower, const std::string cutoffUpper,
    bool wholeImageThreshold);
};


#endif
