///////////////////////////////////////////////////////////////////////////////////////
// fBiasCorrectionDialog.h
//
// Copyright (c) 2020. All rights reserved.
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

#ifndef _fBiasCorrectionDialog_h_
#define _fBiasCorrectionDialog_h_

#include <QDialog>
#include <QObject>

#include "ui_fBiasCorrectionDialog.h"
/**
\class fBiasCorrectionDialog

\brief This class allows users to set parameters for bias correction from the GUI.
*/

class fBiasCorrectionDialog : public QDialog, private Ui::fBiasCorrectionDialog
{
    Q_OBJECT

public:
    fBiasCorrectionDialog();
    ~fBiasCorrectionDialog();

    QString mInputPathName;
    std::string correctionMode;


public slots:
     void CancelButtonPressed();
     void ConfirmButtonPressed();
     void SelectOutputImage();
     void SelectedMode(int);

private:
    void LoadDefaultParameters();

signals:
      void CallBiasCorrection(const std::string correctionType, QString saveFileName,
      int bias_splineOrder, int bias_otsuBins, int bias_maxIterations, int bias_fittingLevels,
      float bias_filterNoise, float bias_fwhm);

};

#endif