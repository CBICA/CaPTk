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

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputT1CEImageButton, SIGNAL(clicked()), this, SLOT(OpenInputT1CEImage()));
  connect(inputT1ImageButton, SIGNAL(clicked()), this, SLOT(OpenInputT1Image()));
  connect(inputT2ImageButton, SIGNAL(clicked()), this, SLOT(OpenInputT2Image()));
  connect(inputFLImageButton, SIGNAL(clicked()), this, SLOT(OpenInputFLImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
}
fBraTSSegmentation::~fBraTSSegmentation()
{
}
void fBraTSSegmentation::CancelButtonPressed()
{
  this->close();
}

void fBraTSSegmentation::OpenInputT1CEImage()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputT1CEImageName->setText(inputImage);

  mInputPathName = cbica::getFilenamePath(inputImage.toStdString(), false).c_str();
}

void fBraTSSegmentation::OpenInputT1Image()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputT1ImageName->setText(inputImage);

  mInputPathName = cbica::getFilenamePath(inputImage.toStdString(), false).c_str();
}

void fBraTSSegmentation::OpenInputT2Image()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputT2ImageName->setText(inputImage);

  mInputPathName = cbica::getFilenamePath(inputImage.toStdString(), false).c_str();
}

void fBraTSSegmentation::OpenInputFLImage()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputFLImageName->setText(inputImage);

  mInputPathName = cbica::getFilenamePath(inputImage.toStdString(), false).c_str();
}

void fBraTSSegmentation::SelectOutputDirectory()
{
  QString outputImage = getExistingDirectory(this, mInputPathName);
  auto temp = outputImage.toStdString();
  if (outputImage.isNull() || outputImage.isEmpty())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = cbica::getFilenamePath(outputImage.toStdString(), false).c_str(); // overwrite previous default path with new output
}

void fBraTSSegmentation::ConfirmButtonPressed()
{
  if (inputT1CEImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the T1CE Image.", this);
    return;
  }
  if (inputT1ImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the T1 Image.", this);
    return;
  }
  if (inputT2ImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the T2 Image.", this);
    return;
  }
  if (inputFLImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the FL Image.", this);
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output file.", this);
    return;
  }

  emit RunBraTSPipeline(
    inputT1CEImageName->text().toStdString(),
    inputT1ImageName->text().toStdString(),
    inputT2ImageName->text().toStdString(),
    inputFLImageName->text().toStdString(),
    outputImageName->text().toStdString());

  this->close();
}
