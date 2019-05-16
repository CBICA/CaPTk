#include "fTexturePipelineDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fTexturePipelineDialog::fTexturePipelineDialog()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
}
fTexturePipelineDialog::~fTexturePipelineDialog()
{
}

void fTexturePipelineDialog::CancelButtonPressed()
{
  this->close();
}

void fTexturePipelineDialog::ConfirmButtonPressed()
{
  auto outputDirName_string = outputDirName->text().toStdString();

  if (outputDirName_string.empty())
  {
    ShowErrorMessage("Please specify the output directory.", this);
    return;
  }

  emit RunTextureFeaturePipeline(outputDirName_string);

  this->close();
}

void fTexturePipelineDialog::SelectOutputDirectory()
{
  QString outputImage = getExistingDirectory(this, mInputPathName);
  if (outputImage.isNull())
    return;
  else
    outputDirName->setText(outputImage);

  QFileInfo fileInfo(outputImage);
  mInputPathName = fileInfo.absoluteFilePath();
}
