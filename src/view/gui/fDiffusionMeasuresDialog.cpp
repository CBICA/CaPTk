#include "fDiffusionMeasuresDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fDiffusionEstimator::fDiffusionEstimator()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(inputMaskButton, SIGNAL(clicked()), this, SLOT(OpenInputMaskImage()));
  connect(inputBvalButton, SIGNAL(clicked()), this, SLOT(OpenInputBValImage()));
  connect(inputBvecButton, SIGNAL(clicked()), this, SLOT(OpenInputBVecImage()));
  connect(inputRegistrationButton, SIGNAL(clicked()), this, SLOT(OpenInputRegistrationImage()));

  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));

  m_ax->setChecked(true);
  m_fa->setChecked(true);
  m_rad->setChecked(true);
  m_tr->setChecked(true);
  m_bzero->setChecked(true);
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
    ShowErrorMessage("Please specify the input image.", this);
    return;
  }
  if ((inputMaskName->text().isEmpty()) || !cbica::isFile(inputMaskName->text().toStdString()))
  {
    ShowErrorMessage("Please specify the mask image.", this);
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
  if (m_register->isChecked() && (inputRegistrationFile->text().isEmpty() || !cbica::isFile(inputRegistrationFile->text().toStdString())))
  {
      ShowErrorMessage("In order to register output, please specify a fixed image.", this);
      return;
  }
	if (outputImageName->text().isEmpty())
	{
    ShowErrorMessage("Please specify the output directory.", this);
		return;
	}
  if (!m_ax->isChecked() && !m_fa->isChecked() && !m_rad->isChecked() && !m_tr->isChecked() && !m_bzero->isChecked())
  {
    ShowErrorMessage("Please select at least one of the given options: Axial diffusivity, fractional anisotropy, radial diffusivity, apparent diffusion coefficient, extract b0");
    return;
  }
  emit RunDiffusionMeasuresCalculation(mInputPathName.toStdString(), mInputMaskName.toStdString(),
      inputBvalName->text().toStdString(), inputBvecName->text().toStdString(),
      m_ax->isChecked(), m_fa->isChecked(), m_rad->isChecked(), m_tr->isChecked(), m_bzero->isChecked(),
    mOutputPathName.toStdString(), m_register->isChecked(), mInputRegistrationFileName.toStdString());

	this->close();
}



void fDiffusionEstimator::OpenInputImage()
{
	auto inputImage = getExistingFile(this, mInputPathName);
	if (inputImage.isNull() || inputImage.isEmpty())
		return;
	else
		inputImageName->setText(inputImage);

	mInputPathName = inputImage;
}

void fDiffusionEstimator::OpenInputMaskImage()
{
	auto inputImage = getExistingFile(this, mInputMaskName);
	if (inputImage.isNull() || inputImage.isEmpty())
		return;
	else
		inputMaskName->setText(inputImage);

	mInputMaskName = inputImage;
}
void fDiffusionEstimator::OpenInputBValImage()
{
  auto inputImage = getExistingFile(this, mInputBValName, "*.bval");
	if (inputImage.isNull() || inputImage.isEmpty())
		return;
	else
		inputBvalName->setText(inputImage);

	mInputBValName = inputImage;
}

void fDiffusionEstimator::OpenInputBVecImage()
{
  auto inputImage = getExistingFile(this, mInputBVecName, "*.bvec");
	if (inputImage.isNull() || inputImage.isEmpty())
		return;
	else
		inputBvecName->setText(inputImage);
	mInputBVecName = inputImage;
}

void fDiffusionEstimator::OpenInputRegistrationImage()
{
    auto inputImage = getExistingFile(this, mInputRegistrationFileName, "*.nii.gz");
    if (inputImage.isNull() || inputImage.isEmpty())
        return;
    else
        inputRegistrationFile->setText(inputImage);
    mInputRegistrationFileName = inputImage;
}

void fDiffusionEstimator::SelectOutputImage()
{
	QString directory = getExistingDirectory(this, mOutputPathName);
	if (directory.isNull() || directory.isEmpty())
		return;
	else
		outputImageName->setText(directory);

	mOutputPathName = directory;
}
