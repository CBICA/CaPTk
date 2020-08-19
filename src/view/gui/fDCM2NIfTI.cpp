#include "fDCM2NIfTI.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fDCM2NIfTIConverter::fDCM2NIfTIConverter()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputDirButton, SIGNAL(clicked()), this, SLOT(OpenInputDir()));
  connect(outputDirButton, SIGNAL(clicked()), this, SLOT(SelectOutputDir()));

  //get OS specific path to dcm2niix binary
  m_exe = getApplicationPath("dcm2niix").c_str();

}
fDCM2NIfTIConverter::~fDCM2NIfTIConverter()
{
}

void fDCM2NIfTIConverter::CancelButtonPressed()
{
  this->close();
}

void fDCM2NIfTIConverter::ConfirmButtonPressed()
{
  auto inputDirName_string = inputDirName->text().toStdString();
  auto outputDirName_string = outputDirName->text().toStdString();

  if ((inputDirName->text().isEmpty()) || !cbica::isDir(inputDirName_string))
  {
    ShowErrorMessage("Please specify the input directory.", this);
    return;
  }
  if (outputDirName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output directory.", this);
    return;
  }

  emit RunDICOMConverter(inputDirName_string, outputDirName_string);
  
  this->close();
}

void fDCM2NIfTIConverter::OpenInputDir()
{
  //auto inputImage = getExistingFile(this, mInputPathName, "DICOM Images (*.dcm *.dicom *.ima *.IMA)");
  auto inputDir = getExistingDirectory(this, mInputPathName);

  if (inputDir.isNull() || inputDir.isEmpty())
    return;
  else
    inputDirName->setText(inputDir);

  mInputPathName = inputDir;
}

void fDCM2NIfTIConverter::SelectOutputDir()
{
  //QString outputImage = getSaveFile(this, mInputPathName, mInputPathName + "dcm2nii.nii.gz");
  QString dir = getExistingDirectory(this, mInputPathName);
  if (dir.isNull() || dir.isEmpty())
    return;
  else
    outputDirName->setText(dir);

  mInputPathName = dir;
}
