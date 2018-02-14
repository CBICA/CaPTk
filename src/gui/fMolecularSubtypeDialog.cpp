#include "fMolecularSubtypeDialog.h"
#include "fProgressDialog.h"

fMolecularSubtypePredictor::fMolecularSubtypePredictor()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
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
fMolecularSubtypePredictor::~fMolecularSubtypePredictor()
{
}
void fMolecularSubtypePredictor::CancelButtonPressed()
{
  this->close();
}
void fMolecularSubtypePredictor::ConfirmButtonPressed()
{
  if (rdExistingClassification->isChecked())
  {
    svmModelButton->setEnabled(true);
    svmModelFileName->setEnabled(true);
    testSubjectsDirectoryButton->setEnabled(true);
    testSubjectsDirectoryName->setEnabled(true);
 
	//--------------SVM model directory--------------------
	if (svmModelFileName->text().toStdString().empty())
	{
		ShowErrorMessage("Please specify the directory of SVM model.");
		return;
	}
	if (!cbica::directoryExists(svmModelFileName->text().toStdString()))
	{
		ShowErrorMessage("The specified SVM model directory does not exist.");
		return;
	}
	//--------------test subjects directory--------------------
	if (testSubjectsDirectoryName->text().toStdString().empty())
	{
		ShowErrorMessage("Please specify the directory of test subjects.");
		return;
	}
	if (!cbica::directoryExists(testSubjectsDirectoryName->text().toStdString()))
	{
		ShowErrorMessage("The specified directory of test subjects does not exist.");
		return;
	}
	//--------------output directory--------------------
	if (outputDirectoryName->text().toStdString().empty())
	{
		ShowErrorMessage("Please specify the output directory.");
		return;
	}
	if (!cbica::directoryExists(outputDirectoryName->text().toStdString()))
	{
		if (!cbica::createDirectory(outputDirectoryName->text().toStdString()))
		{
			ShowErrorMessage("Unable to create the output directory.");
			return;
		}
	}
	//--------------function call--------------------
    emit SurvivalPredictionOnExistingModel(svmModelFileName->text().toStdString(), testSubjectsDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
  else 
  {
    existingMaskDirectoryName->setEnabled(true);
    existingMasksButton->setEnabled(true);

	//--------------input subjects directory--------------------
	if (existingMaskDirectoryName->text().toStdString().empty())
	{
		ShowErrorMessage("Please specify the directory of training subjects.");
		return;
	}
	if (!cbica::directoryExists(existingMaskDirectoryName->text().toStdString()))
	{
		ShowErrorMessage("The specified directory of training subjects does not exist.");
		return;
	}
	//--------------output directory--------------------
	if (outputDirectoryName->text().toStdString().empty())
	{
		ShowErrorMessage("Please specify the output directory.");
		return;
	}
	if (!cbica::directoryExists(outputDirectoryName->text().toStdString()))
	{
		if (!cbica::createDirectory(outputDirectoryName->text().toStdString()))
		{
			ShowErrorMessage("Unable to create the output directory.");
			return;
		}
	}
	//--------------function call--------------------
    emit PrepareNewSurvivalPredictionModel(existingMaskDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
  this->close();
}

void fMolecularSubtypePredictor::OpenSVMModelFile()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    svmModelFileName->setText(directory);

  mInputPathName = directory;
}

void fMolecularSubtypePredictor::OpenExistingMasksDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    existingMaskDirectoryName->setText(directory);

  mInputPathName = directory;
}


void fMolecularSubtypePredictor::SelectOutputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);
}
void fMolecularSubtypePredictor::OpenTestSubjectsDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    testSubjectsDirectoryName->setText(directory);

  mInputPathName = directory;
}

void fMolecularSubtypePredictor::ExistingClassificationRadioButtonChecked()
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
void fMolecularSubtypePredictor::NewModelRadioButtonChecked()
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