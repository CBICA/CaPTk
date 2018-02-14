#include "fDCM2NIfTI.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fDCM2NIfTIConverter::fDCM2NIfTIConverter()
{
  setupUi(this);
  this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));

#if WIN32
  m_exe = "/dcm2nii.exe";

  if (QFile(QApplication::applicationDirPath() + m_exe).exists())
  {
    m_exe = QApplication::applicationDirPath() + m_exe; // populate full application name here itself
  }
  else
  {
    m_exe = QApplication::applicationDirPath() + "/../../src/applications/individualApps/dcm2nii" + m_exe;
  }

#elif __linux__
  m_exe =  "." + QApplication::applicationDirPath() + "/dcm2nii";
#endif

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
  auto inputImageName_string = inputImageName->text().toStdString();
  auto outputImageName_string = outputImageName->text().toStdString();

  if ((inputImageName->text().isEmpty()) || !cbica::isFile(inputImageName_string))
  {
    ShowErrorMessage("Please specify the input Image.", this);
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output file.", this);
    return;
  }

  emit RunDICOMConverter(inputImageName_string, outputImageName_string);
  
  this->close();
}

void fDCM2NIfTIConverter::OpenInputImage()
{
  auto inputImage = getExistingFile(this, mInputPathName, "DICOM Images (*.dcm *.dicom)");

  if (inputImage.isNull())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = inputImage;
}

void fDCM2NIfTIConverter::SelectOutputImage()
{
  QString outputImage = getSaveFile(this, mInputPathName, mInputPathName + "dcm2nii.nii.gz");
  if (outputImage.isNull())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = outputImage;
}
