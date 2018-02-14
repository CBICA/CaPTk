#include "fSkullStripDialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fSkullStripper::fSkullStripper()
{
  setupUi(this);
  this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(referenceImageButton, SIGNAL(clicked()), this, SLOT(OpenReferenceImage()));
  connect(referenceMaskButton, SIGNAL(clicked()), this, SLOT(OpenReferenceMaskImage()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
}
fSkullStripper::~fSkullStripper()
{
}
void fSkullStripper::CancelButtonPressed()
{
  this->close();
}
void fSkullStripper::ConfirmButtonPressed()
{
  auto referenceImageName_string = referenceImageName->text().toStdString();
  auto referenceMaskName_string = referenceMaskName->text().toStdString();
  auto inputImageName_string = inputImageName->text().toStdString();
  auto outputImageName_string = outputImageName->text().toStdString();
  if ((referenceImageName->text().isEmpty()) || !cbica::isFile(referenceImageName_string))
  {
    ShowErrorMessage("Please specify the reference Atlas.");
    return;
  }
  if ((referenceMaskName->text().isEmpty()) || !cbica::isFile(referenceMaskName_string))
  {
    ShowErrorMessage("Please specify the reference Atlas Mask.");
    return;
  }
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

  emit RunSkullStripping(referenceImageName_string, referenceMaskName_string, inputImageName_string, outputImageName_string);
  
  this->close();
}

void fSkullStripper::OpenReferenceImage()
{
  auto referenceAtlas = getExistingFile(this, mInputPathName);
  if (referenceAtlas.isNull())
    return;
  else
    referenceImageName->setText(referenceAtlas);

  mInputPathName = referenceAtlas;
}

void fSkullStripper::OpenReferenceMaskImage()
{
  auto referenceMask = getExistingFile(this, mInputPathName);
  if (referenceMask.isNull())
    return;
  else
    referenceMaskName->setText(referenceMask);

  mInputPathName = referenceMask;
}

void fSkullStripper::OpenInputImage()
{
  auto inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = inputImage;
}

void fSkullStripper::SelectOutputImage()
{
  QString outputImage;
  if (inputImageName->text().isEmpty())
  {
    outputImage = getSaveFile(this, mInputPathName);
  }
  else
  {
    outputImage = getSaveFile(this, mInputPathName, mInputPathName + "skullStrip_Output.nii.gz");
  }

  if (outputImage.isNull())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = outputImage;
}
