///////////////////////////////////////////////////////////////////////////////////////
// fDiffusionEstimator.h
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

#ifndef _fDiffusionMeasuresDialog_h_
#define _fDiffusionMeasuresDialog_h_


//#include "CAPTk.h"
#include "ui_fDiffusionMeasuresDialog.h"

/**
\class fDiffusionEstimator

\brief This class controls the elements in the recurrence dialog
*/
class fDiffusionEstimator : public QDialog, private Ui::fDiffusionEstimator
{
	Q_OBJECT

public:
	fDiffusionEstimator();
	~fDiffusionEstimator();
	int mode;

	void SetCurrentImagePath(const QString &inputPath)
	{
		mInputPathName = inputPath;
    outputImageName->setText(mInputPathName);
	}

	QString mInputPathName;
	QString mInputMaskName;
	QString mInputBValName;
	QString mInputBVecName;
	QString mOutputPathName;


	public slots:
	void CancelButtonPressed();
	void ConfirmButtonPressed();
	void OpenInputImage();
	void OpenInputMaskImage();
	void OpenInputBValImage();
	void OpenInputBVecImage();
	void SelectOutputImage();

signals:
	void RunDiffusionMeasuresCalculation(const std::string inputImageFile, const std::string inputMaskFile,
    const std::string bValFile, const std::string bVecFile, const bool ax, const bool fa, const bool rad, const bool tr, const std::string outputFolder);
};


#endif
