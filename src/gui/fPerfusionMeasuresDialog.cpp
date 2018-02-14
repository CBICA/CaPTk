#include "fPerfusionMeasuresDialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fPerfusionEstimator::fPerfusionEstimator()
{
  setupUi(this);
  this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));

  m_rcbv->setChecked(true);
  m_psr->setChecked(true);
  m_ph->setChecked(true);
}
fPerfusionEstimator::~fPerfusionEstimator()
{
}
void fPerfusionEstimator::CancelButtonPressed()
{
  this->close();
}
void fPerfusionEstimator::ConfirmButtonPressed()
{
  auto inputImageName_string = inputImageName->text().toStdString();
  auto outputImageName_string = outputImageName->text().toStdString();

  if ((inputImageName->text().isEmpty()) || !cbica::isFile(inputImageName_string))
  {
    ShowErrorMessage("Please specify the input Image.");
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output file.");
    return;
  }
  emit RunPerfusionMeasuresCalculation(mInputPathName.toStdString(), m_rcbv->isChecked(),m_psr->isChecked(), m_ph->isChecked(), mOutputPathName.toStdString());
  
  this->close();
}



void fPerfusionEstimator::OpenInputImage()
{
  auto inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = inputImage;
}

void fPerfusionEstimator::SelectOutputImage()
{
	QString directory = getExistingDirectory(this, mOutputPathName);
	if (directory.isNull())
		return;
	else
		outputImageName->setText(directory);

	mOutputPathName = directory;
}
