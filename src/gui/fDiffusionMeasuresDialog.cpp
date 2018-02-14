#include "fDiffusionMeasuresDialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fDiffusionEstimator::fDiffusionEstimator()
{
  setupUi(this);
  this->setModal(true); // this is a pre-processing routine and therefore should be modal
	this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

	connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
	connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
	connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
	connect(inputMaskButton, SIGNAL(clicked()), this, SLOT(OpenInputMaskImage()));
	connect(inputBvalButton, SIGNAL(clicked()), this, SLOT(OpenInputBValImage()));
	connect(inputBvecButton, SIGNAL(clicked()), this, SLOT(OpenInputBVecImage()));


	connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));

	m_ax->setChecked(true);
	m_fa->setChecked(true);
	m_rad->setChecked(true);
	m_tr->setChecked(true);
}
fDiffusionEstimator::~fDiffusionEstimator()
{
}
void fDiffusionEstimator::CancelButtonPressed()
{
	this->close();
}
void fDiffusionEstimator::ConfirmButtonPressed()
{
	auto inputImageName_string = inputImageName->text().toStdString();
	auto outputImageName_string = outputImageName->text().toStdString();

	if ((inputImageName->text().isEmpty()) || !cbica::isFile(inputImageName_string))
	{
    ShowErrorMessage("Please specify the input Image.", this);
		return;
	}
	if ((inputMaskName->text().isEmpty()) || !cbica::isFile(inputMaskName->text().toStdString()))
	{
    ShowErrorMessage("Please specify the mask Image.", this);
		return;
	}
	if ((inputBvalName->text().isEmpty()) || !cbica::isFile(inputBvalName->text().toStdString()))
	{
    ShowErrorMessage("Please specify the BVal file.", this);
		return;
	}
	if ((inputBvecName->text().isEmpty()) || !cbica::isFile(inputBvecName->text().toStdString()))
	{
    ShowErrorMessage("Please specify the BVec file.", this);
		return;
	}
	if (outputImageName->text().isEmpty())
	{
    ShowErrorMessage("Please specify the output file.", this);
		return;
	}
	emit RunDiffusionMeasuresCalculation(mInputPathName.toStdString(), mInputMaskName.toStdString(),mInputBValName.toStdString(),mInputBVecName.toStdString(),m_ax->isChecked(), m_fa->isChecked(), m_rad->isChecked(), m_tr->isChecked(), mOutputPathName.toStdString());

	this->close();
}



void fDiffusionEstimator::OpenInputImage()
{
	auto inputImage = getExistingFile(this, mInputPathName);
	if (inputImage.isNull())
		return;
	else
		inputImageName->setText(inputImage);

	mInputPathName = inputImage;
}

void fDiffusionEstimator::OpenInputMaskImage()
{
	auto inputImage = getExistingFile(this, mInputMaskName);
	if (inputImage.isNull())
		return;
	else
		inputMaskName->setText(inputImage);

	mInputMaskName = inputImage;
}
void fDiffusionEstimator::OpenInputBValImage()
{
	auto inputImage = getExistingFile(this, mInputBValName);
	if (inputImage.isNull())
		return;
	else
		inputBvalName->setText(inputImage);

	mInputBValName = inputImage;
}

void fDiffusionEstimator::OpenInputBVecImage()
{
	auto inputImage = getExistingFile(this, mInputBVecName);
	if (inputImage.isNull())
		return;
	else
		inputBvecName->setText(inputImage);
	mInputBVecName = inputImage;
}

void fDiffusionEstimator::SelectOutputImage()
{
	QString directory = getExistingDirectory(this, mOutputPathName);
	if (directory.isNull())
		return;
	else
		outputImageName->setText(directory);

	mOutputPathName = directory;
}
