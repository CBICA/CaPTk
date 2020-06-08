#include "fBraTSSegmentation.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fBraTSSegmentation::fBraTSSegmentation()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(referenceImageButton, SIGNAL(clicked()), this, SLOT(OpenReferenceImage()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
}
fBraTSSegmentation::~fBraTSSegmentation()
{
}
void fBraTSSegmentation::CancelButtonPressed()
{
  this->close();
}
void fBraTSSegmentation::ConfirmButtonPressed()
{
  if (referenceImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the reference image.", this);
    return;
  }
  if (inputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the input Image.", this);
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output file.", this);
    return;
  }

  emit RunHistogramMatching(referenceImageName->text().toStdString(), inputImageName->text().toStdString(), outputImageName->text().toStdString());
  
  this->close();
}

void fBraTSSegmentation::OpenReferenceImage()
{
  QString referenceAtlas = getExistingFile(this, mInputPathName);
  if (referenceAtlas.isNull() || referenceAtlas.isEmpty())
    return;
  else
    referenceImageName->setText(referenceAtlas);

  mInputPathName = cbica::getFilenamePath(referenceAtlas.toStdString(), false).c_str();
}

void fBraTSSegmentation::OpenInputImage()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = cbica::getFilenamePath(inputImage.toStdString(), false).c_str();
}

void fBraTSSegmentation::SelectOutputImage()
{
  QString outputImage = getSaveFile(this, mInputPathName, mInputPathName + "histoMatch_output.nii.gz");
  auto temp = outputImage.toStdString();
  if (outputImage.isNull() || outputImage.isEmpty())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = cbica::getFilenamePath(outputImage.toStdString(), false).c_str(); // overwrite previous default path with new output
}
