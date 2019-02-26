#include "fDeepMedicNormDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fDeepMedicNormalizer::fDeepMedicNormalizer()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(maskImageButton, SIGNAL(clicked()), this, SLOT(OpenMaskImage()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
}
fDeepMedicNormalizer::~fDeepMedicNormalizer()
{
}
void fDeepMedicNormalizer::CancelButtonPressed()
{
  this->close();
}
void fDeepMedicNormalizer::ConfirmButtonPressed()
{
  if (maskImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the mask image.", this);
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

  emit RunDeepMedicNormalizer(inputImageName->text().toStdString(), maskImageName->text().toStdString(), outputImageName->text().toStdString(),
    options_quantileLowerName->text().toStdString(), options_quantileUpperName->text().toStdString(),
    options_cutoffLowerName->text().toStdString(), options_cutoffUpperName->text().toStdString(),
    wholeImageThresholdCheckBoxBox->isChecked());
  
  this->close();
}

void fDeepMedicNormalizer::OpenMaskImage()
{
  QString maskImage = getExistingFile(this, mInputPathName);
  if (maskImage.isNull() || maskImage.isEmpty())
    return;
  else
    maskImageName->setText(maskImage);

  mInputPathName = cbica::getFilenameBase(maskImage.toStdString(), false).c_str();
}

void fDeepMedicNormalizer::OpenInputImage()
{
  QString inputImage = getExistingFile(this, mInputPathName);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = cbica::getFilenameBase(inputImage.toStdString(), false).c_str();
}

void fDeepMedicNormalizer::SelectOutputImage()
{
  QString outputImage = getSaveFile(this, mInputPathName, outputImageName->text());
  if (outputImage.isNull() || outputImage.isEmpty())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = cbica::getFilenameBase(outputImage.toStdString(), false).c_str(); // overwrite previous default path with new output
}
