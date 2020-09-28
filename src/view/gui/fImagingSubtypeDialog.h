///////////////////////////////////////////////////////////////////////////////////////
// fImagingSubtypePredictor.h
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

#ifndef _fImagingSubtypePredictor_h_
#define _fImagingSubtypePredictor_h_


//#include "CAPTk.h"
#include "ui_fImagingSubtypeDialog.h"

#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fImagingSubtypePredictor

\brief This class controls the elements in the recurrence dialog
*/
class fImagingSubtypePredictor : public QDialog, private Ui::fImagingSubtypePredictor
{
	Q_OBJECT

public:
	fImagingSubtypePredictor();
	~fImagingSubtypePredictor();
	int mode;

	void SetCurrentImagePath(const QString &inputPath)
	{
		mInputPathName = inputPath;
	}

	QString mInputPathName;


	public slots:
	void CancelButtonPressed();
	void ConfirmButtonPressed();
	void OpenSVMModelFile();
	void OpenExistingMasksDirectory();
	void SelectOutputDirectory();
	void ExistingClassificationRadioButtonChecked();
	void NewModelRadioButtonChecked();
	void OpenTestSubjectsDirectory();

signals:
	void PrepareNewSurvivalPredictionModel(const std::string inputdirectory, const std::string outputdirectory);
	void SurvivalPredictionOnExistingModel(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);
	void SubjectBasedSurvivalEstimate(const std::string output, const std::string model, double age);
};


#endif
