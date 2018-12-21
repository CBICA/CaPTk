#include "fPreprocessingDialog.h"
#include "fProgressDialog.h"

#include "qobject.h"

fPreprocessingDialog::fPreprocessingDialog()
{
  setupUi(this);
  mode = -1;
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(existingMasksButton, SIGNAL(clicked()), this, SLOT(OpenExistingMasksDirectory()));
  connect(svmModelButton, SIGNAL(clicked()), this, SLOT(OpenSVMModelFile()));
  connect(testSubjectsDirectoryButton, SIGNAL(clicked()), this, SLOT(OpenTestSubjectsDirectory()));
  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
  connect(rdExistingClassification, SIGNAL(toggled(bool)), this, SLOT(ExistingClassificationRadioButtonChecked()));
  connect(rdCreateModel, SIGNAL(toggled(bool)), this, SLOT(NewModelRadioButtonChecked()));

  rdNewClassification->setChecked(true);
  rdExistingClassification->setChecked(false);
  svmModelButton->setEnabled(false);
  svmModelFileName->setEnabled(false);
  testSubjectsDirectoryButton->setEnabled(false);
  testSubjectsDirectoryName->setEnabled(false);

  existingMaskDirectoryName->setEnabled(false);
  existingMasksButton->setEnabled(false);
}
fPreprocessingDialog::~fPreprocessingDialog()
{
}

void fPreprocessingDialog::CancelButtonPressed()
{
  this->close();
}

void fPreprocessingDialog::ConfirmButtonPressed()
{

}

void fPreprocessingDialog::OpenExistingMasksDirectory()
{

}

void fPreprocessingDialog::OpenSVMModelFile()
{

}

void fPreprocessingDialog::OpenTestSubjectsDirectory()
{

}

void fPreprocessingDialog::SelectOutputDirectory()
{

}

void fPreprocessingDialog::ExistingClassificationRadioButtonChecked()
{

}

void fPreprocessingDialog::NewModelRadioButtonChecked()
{

}
