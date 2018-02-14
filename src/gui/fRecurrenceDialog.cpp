#include "fRecurrenceDialog.h"
#include "fProgressDialog.h"

fRecurrenceDialog::fRecurrenceDialog()
{
  this->setWindowTitle("Glioblastoma Infiltration Index");
  //TBD read this from common place 
  setupUi(this);
  //this->setWindowModality(Qt::NonModal);
  this->setModal(true);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  //this->setFixedWidth(400);
  //this->setFixedHeight(300);
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(existingMasksButton, SIGNAL(clicked()), this, SLOT(OpenExistingMasksDirectory()));
  connect(svmModelButton, SIGNAL(clicked()), this, SLOT(OpenSVMModelFile()));
  connect(testSubjectsDirectoryButton, SIGNAL(clicked()), this, SLOT(OpenTestSubjectsDirectory()));
  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
  connect(rdExistingClassification, SIGNAL(toggled(bool)), this, SLOT(ExistingClassificationRadioButtonChecked()));
  connect(rdCreateModel, SIGNAL(toggled(bool)), this, SLOT(NewModelRadioButtonChecked()));
  connect(rdNewClassification, SIGNAL(toggled(bool)), this, SLOT(CurrentSubjectRadioButtonChecked()));

  rdNewClassification->setChecked(true);
  rdExistingClassification->setChecked(false);
  svmModelButton->setEnabled(false);
  svmModelFileName->setEnabled(false);
  testSubjectsDirectoryButton->setEnabled(false);
  testSubjectsDirectoryName->setEnabled(false);
  existingMaskDirectoryName->setEnabled(false);
  existingMasksButton->setEnabled(false);
  cbT1Data->hide();
  cbDistanceData->hide();
  cbDTIData->hide();
  cbPerfData->hide();
}
fRecurrenceDialog::~fRecurrenceDialog()
{
}
void fRecurrenceDialog::CancelButtonPressed()
{
  this->close();
}
void fRecurrenceDialog::ConfirmButtonPressed()
{
	if (outputDirectoryName->text().isEmpty())
	{
		ShowErrorMessage("Please specify the output directory.");
		return;
	}
	if (!cbica::directoryExists(outputDirectoryName->text().toStdString()))
	{
		if (!cbica::createDirectory(outputDirectoryName->text().toStdString()))
		{
			ShowErrorMessage("The output directory can not be created.");
			return;
		}
	}

  if (rdNewClassification->isChecked())
  {
   // emit SubjectBasedRecurrenceEstimate(outputDirectoryName->text().toStdString(), cbT1Data->isChecked(), cbDTIData->isChecked(), cbPerfData->isChecked(), cbDistanceData->isChecked());
	  emit SubjectBasedRecurrenceEstimate(outputDirectoryName->text().toStdString(), true,true,true,true);
    this->close();
  }
  else if (rdLoadedClassification->isChecked())
  {
	  // emit SubjectBasedRecurrenceEstimate(outputDirectoryName->text().toStdString(), cbT1Data->isChecked(), cbDTIData->isChecked(), cbPerfData->isChecked(), cbDistanceData->isChecked());
	  emit SubjectBasedExistingRecurrenceEstimate(outputDirectoryName->text().toStdString(), true, true, true, true);
	  this->close();
  }
  else if (rdExistingClassification->isChecked())
  {
    if (svmModelFileName->text().isEmpty())
    {
      ShowErrorMessage("Please specify the directory of SVM model.");
      return;
    }
    if (testSubjectsDirectoryName->text().isEmpty())
    {
      ShowErrorMessage("Please specify the directory of test subjects.");
      return;
    }

    //if (!cbT1Data->isChecked() && !cbDTIData->isChecked() && !cbPerfData->isChecked() && !cbDistanceData->isChecked())
    //{
    //  ShowErrorMessage("Please specify the required modalities.");
    //  return;
    //}
    //emit ExistingModelBasedRecurrenceEstimate(svmModelFileName->text().toStdString(), testSubjectsDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString(), cbT1Data->isChecked(), cbDTIData->isChecked(), cbPerfData->isChecked(), cbDistanceData->isChecked());
	emit ExistingModelBasedRecurrenceEstimate(svmModelFileName->text().toStdString(), testSubjectsDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString(), true,true,true,true);
  }
  else
    //emit TrainNewModel(existingMaskDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString(), cbT1Data->isChecked(), cbDTIData->isChecked(), cbPerfData->isChecked(), cbDistanceData->isChecked());
	emit TrainNewModel(existingMaskDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString(), true,true,true,true);
  this->close();
}

void fRecurrenceDialog::OpenSVMModelFile()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    svmModelFileName->setText(directory);

  mInputPathName = directory;
}

void fRecurrenceDialog::OpenExistingMasksDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    existingMaskDirectoryName->setText(directory);

  mInputPathName = directory;
}


void fRecurrenceDialog::SelectOutputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);
}
void fRecurrenceDialog::OpenTestSubjectsDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    testSubjectsDirectoryName->setText(directory);

  mInputPathName = directory;
}

void fRecurrenceDialog::ExistingClassificationRadioButtonChecked()
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
void fRecurrenceDialog::NewModelRadioButtonChecked()
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
void fRecurrenceDialog::CurrentSubjectRadioButtonChecked()
{
	if (rdNewClassification->isChecked())
	{
		svmModelButton->setEnabled(false);
		svmModelFileName->setEnabled(false);
		testSubjectsDirectoryButton->setEnabled(false);
		testSubjectsDirectoryName->setEnabled(false);
		existingMaskDirectoryName->setEnabled(false);
		existingMasksButton->setEnabled(false);
	}
}