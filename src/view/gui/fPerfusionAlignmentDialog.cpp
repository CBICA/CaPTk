#include "fPerfusionAlignmentDialog.h"
#include "CaPTkGUIUtils.h"
#include "cbicaITKUtilities.h"

fPerfusionAligner::fPerfusionAligner()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(SelectInputImage()));
}
fPerfusionAligner::~fPerfusionAligner()
{
}
void fPerfusionAligner::CancelButtonPressed()
{
  this->close();
}
void fPerfusionAligner::ConfirmButtonPressed()
{
  if ((inputImageName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the DSC-MRI Image.");
    return;
  }
  if ((inputBeforePointsName->text().isEmpty()))
  {
	  ShowErrorMessage("Please specify the number of points to pick before the drop.");
	  return;
  }
  if ((inputAfterPointsName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the number of points to pick after the drop.");
    return;
  }
  if ((inputEchoTimeName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the echo time.");
    return;
  }
  if ((inputBaselineLine->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the echo time.");
    return;
  }
  if ((inputScaleDropBeforeMeanName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the echo time.");
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output folder.");
    return;
  }
  emit RunPerfusionAlignmentCalculation(inputEchoTimeName->text().toDouble(), outputEchoTimeName->text().toDouble(), inputBeforePointsName->text().toInt(), inputAfterPointsName->text().toInt(), inputScaleDropBeforeMeanName->text().toInt(),
    inputBaselineLine->text().toInt(), inputImageName->text().toStdString(), inputMaskName->text().toStdString(), mOutputPathName.toStdString());

  this->close();
}

void fPerfusionAligner::SelectOutputImage()
{
  QString directory = getExistingDirectory(this, mOutputPathName);
  if (directory.isNull() || directory.isEmpty())
    return;
  else
    outputImageName->setText(directory);

  mOutputPathName = directory;
}

void fPerfusionAligner::SelectInputImage()
{
	auto inputImage = getExistingFile(this, mInputPathName);
	if (inputImage.isNull() || inputImage.isEmpty())
		return;
	else
		inputImageName->setText(inputImage);

	mInputPathName = inputImage;
}