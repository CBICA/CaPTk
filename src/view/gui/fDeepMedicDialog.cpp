#include "fDeepMedicDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fDeepMedicDialog::fDeepMedicDialog()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
  connect(modelImageButton, SIGNAL(clicked()), this, SLOT(SelectModelDirectory()));
  connect(brainTumorSegmentationButton, SIGNAL(toggled(bool)), this, SLOT(SetDefaultModel()));
  connect(skullStrippingButton, SIGNAL(toggled(bool)), this, SLOT(SetDefaultModel()));
  connect(customButton, SIGNAL(toggled(bool)), this, SLOT(SetDefaultModel()));
  //connect(brainTumorSegmentationButton, SIGNAL(toggled(bool)), this, [this] { SetDefaultModel(fDeepMedicDialog::Tumor); });

  outputDirName->setText(mInputPathName);
  m_baseModelDir = cbica::normPath(getCaPTkDataDir() + "/deepMedic/saved_models/");

  m_exe = getApplicationPath("DeepMedic").c_str();

  m_dataDir = getCaPTkDataDir().c_str();
}
fDeepMedicDialog::~fDeepMedicDialog()
{
}

void fDeepMedicDialog::CancelButtonPressed()
{
  this->close();
}

void fDeepMedicDialog::SetDefaultModel()
{
  if (brainTumorSegmentationButton->isChecked())
  {
    SetDefaultModel(fDeepMedicDialog::Tumor);
  }
  else if (skullStrippingButton->isChecked())
  {
    SetDefaultModel(fDeepMedicDialog::SkullStripping);
  }
  else if (customButton->isChecked())
  {
    SetDefaultModel(fDeepMedicDialog::Custom);
  }
}

void fDeepMedicDialog::SetDefaultModel(int modelType)
{
  if (modelType == fDeepMedicDialog::Tumor)
  {
    brainTumorSegmentationButton->setChecked(true);
    skullStrippingButton->setChecked(false);
    customButton->setChecked(false);
    modelImageButton->setChecked(false);
    auto currentModelDir = m_baseModelDir  + "/brainTumorSegmentation/";
    modelDirName->setReadOnly(true);
    modelDirName->setText(currentModelDir.c_str());
  }
  else if (modelType == fDeepMedicDialog::SkullStripping)
  {
    brainTumorSegmentationButton->setChecked(false);
    skullStrippingButton->setChecked(true);
    customButton->setChecked(false);
    modelImageButton->setChecked(false);
    auto currentModelDir = m_baseModelDir + "/skullStripping/";
    modelDirName->setReadOnly(true);
    modelDirName->setText(currentModelDir.c_str());
  }
  else if (modelType == fDeepMedicDialog::Custom)
  {
    brainTumorSegmentationButton->setChecked(false);
    skullStrippingButton->setChecked(false);
    customButton->setChecked(true);
    modelImageButton->setChecked(true);
    auto currentModelDir = m_baseModelDir;
    modelDirName->setReadOnly(false);
    modelDirName->setText(currentModelDir.c_str());
  }
}

void fDeepMedicDialog::ConfirmButtonPressed()
{
  auto outputDirName_string = outputDirName->text().toStdString();
  auto modelDirName_string = modelDirName->text().toStdString();

  if (outputDirName_string.empty())
  {
    ShowErrorMessage("Please specify the output directory.", this);
    return;
  }
  if (modelDirName_string.empty())
  {
    ShowErrorMessage("Please specify the model directory.", this);
    return;
  }

  emit RunDeepMedic(modelDirName_string, outputDirName_string);

  this->close();
}

void fDeepMedicDialog::SelectOutputDirectory()
{
  QString outputImage = getExistingDirectory(this, mInputPathName);
  if (outputImage.isNull())
    return;
  else
    outputDirName->setText(outputImage);

  QFileInfo fileInfo(outputImage);
  mInputPathName = fileInfo.absoluteFilePath();
}

void fDeepMedicDialog::SelectModelDirectory()
{
  QString modelDir = getExistingDirectory(this, modelDirName->text());
  if (modelDir.isNull())
    return;
  else
    modelDirName->setText(modelDir);
}
