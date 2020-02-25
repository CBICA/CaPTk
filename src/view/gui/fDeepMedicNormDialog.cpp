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
  std::string defaultExtension = ".nii.gz";
  std::string path, base, ext; // reusable variables for split filenames
  std::string outputWithExtension = outputImageName->text().toStdString();
  // Safety check, add a .nii.gz to the default filename
  cbica::splitFileName(outputImageName->text().toStdString(), path, base, ext);
  if (ext.empty()) {
      outputWithExtension += defaultExtension;
  }

  QString outputImage = getSaveFile(this, mInputPathName, QString::fromStdString(outputWithExtension));
  cbica::splitFileName(outputImage.toStdString(), path, base, ext);
  if (!base.empty() && ext.empty()) {
      // if the user deliberately specifies an output file with no extension,
      // it will cause a crash due to how ITK handles file writing.
      // We could raise an error here asking for an extension, but for now just append .nii.gz
      outputImage += QString::fromStdString(defaultExtension);

  }
  if (outputImage.isNull() || outputImage.isEmpty())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = cbica::getFilenamePath(outputImage.toStdString(), false).c_str(); // overwrite previous default path with new output
}
