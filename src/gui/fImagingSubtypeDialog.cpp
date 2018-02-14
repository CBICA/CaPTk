#include "fImagingSubtypeDialog.h"
#include "fProgressDialog.h"

fImagingSubtypePredictor::fImagingSubtypePredictor()
{
  setupUi(this);
  //this->setWindowModality(Qt::NonModal);
  this->setModal(true);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(existingMasksButton, SIGNAL(clicked()), this, SLOT(OpenExistingMasksDirectory()));
  connect(svmModelButton, SIGNAL(clicked()), this, SLOT(OpenSVMModelFile()));
  connect(testSubjectsDirectoryButton, SIGNAL(clicked()), this, SLOT(OpenTestSubjectsDirectory()));
  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
  connect(rdExistingClassification, SIGNAL(toggled(bool)), this, SLOT(ExistingClassificationRadioButtonChecked()));
  connect(rdCreateModel, SIGNAL(toggled(bool)), this, SLOT(NewModelRadioButtonChecked()));

  //rdExistingClassification->setEnabled(false);
  //rdExistingClassification->setChecked(true);
  rdCreateModel->setChecked(true);
  NewModelRadioButtonChecked();
  
  //existingMaskDirectoryName->setEnabled(false);
  //existingMasksButton->setEnabled(false);
}
fImagingSubtypePredictor::~fImagingSubtypePredictor()
{
}
void fImagingSubtypePredictor::CancelButtonPressed()
{
  this->close();
}
void fImagingSubtypePredictor::ConfirmButtonPressed()
{
  if (rdExistingClassification->isChecked())
  {
    svmModelButton->setEnabled(true);
    svmModelFileName->setEnabled(true);
    testSubjectsDirectoryButton->setEnabled(true);
    testSubjectsDirectoryName->setEnabled(true);
    if (svmModelFileName->text().toStdString().empty())
    {
      ShowErrorMessage("Please specify the directory of SVM model.", this);
      return;
    }
    if (testSubjectsDirectoryName->text().toStdString().empty())
    {
      ShowErrorMessage("Please specify the directory of test subjects.", this);
      return;
    }
    if (outputDirectoryName->text().toStdString().empty())
    {
      ShowErrorMessage("Please specify the output directory.", this);
      return;
    }
    emit SurvivalPredictionOnExistingModel(svmModelFileName->text().toStdString(), testSubjectsDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
  else 
  {
    existingMaskDirectoryName->setEnabled(true);
    existingMasksButton->setEnabled(true);
    emit PrepareNewSurvivalPredictionModel(existingMaskDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
  this->close();
}

void fImagingSubtypePredictor::OpenSVMModelFile()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    svmModelFileName->setText(directory);

  mInputPathName = directory;
}

void fImagingSubtypePredictor::OpenExistingMasksDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    existingMaskDirectoryName->setText(directory);

  mInputPathName = directory;
}


void fImagingSubtypePredictor::SelectOutputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);
}
void fImagingSubtypePredictor::OpenTestSubjectsDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    testSubjectsDirectoryName->setText(directory);

  mInputPathName = directory;
}

void fImagingSubtypePredictor::ExistingClassificationRadioButtonChecked()
{
  if (rdExistingClassification->isChecked())
  {
    svmModelButton->setEnabled(true);
    svmModelFileName->setEnabled(true);
    testSubjectsDirectoryButton->setEnabled(true);
    testSubjectsDirectoryName->setEnabled(true);
    existingMaskDirectoryName->setEnabled(false);
    existingMasksButton->setEnabled(false);
  }
}
void fImagingSubtypePredictor::NewModelRadioButtonChecked()
{
  if (rdCreateModel->isChecked())
  {
    svmModelButton->setEnabled(false);
    svmModelFileName->setEnabled(false);
    testSubjectsDirectoryButton->setEnabled(false);
    testSubjectsDirectoryName->setEnabled(false);
    existingMaskDirectoryName->setEnabled(true);
    existingMasksButton->setEnabled(true);
  }
}