#include "fDirectionalityDialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fDirectionalityDialog::fDirectionalityDialog()
{
  setupUi(this);
  //this->setWindowModality(Qt::NonModal);
  this->setModal(true);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(roi1ImageButton, SIGNAL(clicked()), this, SLOT(OpenROI1Image()));
  connect(roi2ImageButton, SIGNAL(clicked()), this, SLOT(OpenROI2Image()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputFolder()));
}
fDirectionalityDialog::~fDirectionalityDialog()
{
}
void fDirectionalityDialog::CancelButtonPressed()
{
  this->close();
}
void fDirectionalityDialog::ConfirmButtonPressed()
{
  if (roi1ImageName->text().toStdString().empty())
  {
    ShowErrorMessage("Please specify the reference image.", this);
    return;
  }
  if (roi1ImageName->text().toStdString().empty())
  {
    ShowErrorMessage("Please specify the input Image.", this);
    return;
  }
  if (outputDirectory->text().toStdString().empty())
  {
    ShowErrorMessage("Please specify the output file.", this);
    return;
  }

  emit RunDirectionalityEstimator(roi1ImageName->text().toStdString(), roi2ImageName->text().toStdString(), outputDirectory->text().toStdString());
  
  this->close();
}

void fDirectionalityDialog::OpenROI1Image()
{
  QString roi1File = getExistingFile(this, mInputPathName);
  if (roi1File.isNull())
    return;
  else
    roi1ImageName->setText(roi1File);

  mInputPathName = cbica::getFilenamePath(roi1File.toStdString(), false).c_str();
}

void fDirectionalityDialog::OpenROI2Image()
{
  QString roi2File = getExistingFile(this, mInputPathName);
  if (roi2File.isNull())
    return;
  else
    roi2ImageName->setText(roi2File);

  mInputPathName = cbica::getFilenamePath(roi2File.toStdString(), false).c_str();
}

void fDirectionalityDialog::SelectOutputFolder()
{
  QString outputFolder = getExistingDirectory(this, mInputPathName);
  if (outputFolder.isNull())
    return;
  else
    outputDirectory->setText(outputFolder);

  mInputPathName = cbica::getFilenamePath(outputFolder.toStdString(), false).c_str(); // overwrite previous default path with new output
}
