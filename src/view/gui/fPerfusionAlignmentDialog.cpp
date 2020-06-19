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
  connect(inputT1ceImageButton, SIGNAL(clicked()), this, SLOT(SelectT1ceInputImage()));
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
  if ((inputT1ceImageName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the T1ce Image.");
    return;
  }
  if ((inputBeforePointsLabel->text().isEmpty()))
  {
	  ShowErrorMessage("Please specify the number of points to pick before the drop.");
	  return;
  }
  if ((inputAfterPointsLabel->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the number of points to pick after the drop.");
    return;
  }
  if ((inputEchoTimeName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the echo time.");
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output folder.");
    return;
  }
  emit RunPerfusionAlignmentCalculation(inputEchoTimeName->text().toInt(),inputBeforePointsName->text().toInt(), inputAfterPointsName->text().toInt(), mInputPathName.toStdString(), mInputT1cePathName.toStdString(), mOutputPathName.toStdString());

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
void fPerfusionAligner::SelectT1ceInputImage()
{
  auto inputT1ceImage = getExistingFile(this, mInputT1cePathName);
  if (inputT1ceImage.isNull() || inputT1ceImage.isEmpty())
    return;
  else
    inputT1ceImageName->setText(inputT1ceImage);
  mInputT1cePathName = inputT1ceImage;
}