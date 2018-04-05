#include "fHistoMatchDialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fHistoMatcher::fHistoMatcher()
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
fHistoMatcher::~fHistoMatcher()
{
}
void fHistoMatcher::CancelButtonPressed()
{
  this->close();
}
void fHistoMatcher::ConfirmButtonPressed()
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

void fHistoMatcher::OpenReferenceImage()
{
  QString referenceAtlas = getExistingFile(this, mInputPathName);
  if (referenceAtlas.isNull())
    return;
  else
    referenceImageName->setText(referenceAtlas);

  mInputPathName = cbica::getFilenameBase(referenceAtlas.toStdString(), false).c_str();
}

void fHistoMatcher::OpenInputImage()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = cbica::getFilenameBase(inputImage.toStdString(), false).c_str();
}

void fHistoMatcher::SelectOutputImage()
{
  QString outputImage = getSaveFile(this, mInputPathName, mInputPathName + "histoMatch_output.nii.gz");
  if (outputImage.isNull())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = cbica::getFilenameBase(outputImage.toStdString(), false).c_str(); // overwrite previous default path with new output
}
