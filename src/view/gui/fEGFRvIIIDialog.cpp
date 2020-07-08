#include "fEGFRvIIIDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "cbicaLogging.h"

fEGFRvIIIPredictor::fEGFRvIIIPredictor()
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
  connect(disclaimerButton, SIGNAL(clicked()), this, SLOT(CheckForDisclaimer()));

  disclaimerLabel->setText("You can find a pretrained model, based on a study at PENN, ");
  //rdExistingClassification->setEnabled(false);
  rdExistingClassification->setChecked(true); // force testing
  rdExistingClassification->setVisible(false); // hide radio button

  ExistingClassificationRadioButtonChecked();
  
  // Explicitly disable training.
  existingMaskDirectoryName->setEnabled(false);
  existingMasksButton->setEnabled(false);
  rdCreateModel->setEnabled(false);
  existingMaskDirectoryName->setVisible(false);
  existingMasksButton->setVisible(false);
  rdCreateModel->setVisible(false);
  trainingDirectoryLabel->setVisible(false);

}
fEGFRvIIIPredictor::~fEGFRvIIIPredictor()
{
}
void fEGFRvIIIPredictor::CancelButtonPressed()
{
  this->close();
}
void fEGFRvIIIPredictor::ConfirmButtonPressed()
{
  if (rdExistingClassification->isChecked())
  {
    svmModelButton->setEnabled(true);
    svmModelFileName->setEnabled(true);
    testSubjectsDirectoryButton->setEnabled(true);
    testSubjectsDirectoryName->setEnabled(true);

	//--------------SVM model directory--------------------
    //if (svmModelFileName->text().isEmpty())
    //{
    //  ShowErrorMessage("Please specify the directory of SVM model.");
    //  return;
    //}
    //if (!cbica::directoryExists(svmModelFileName->text().toStdString()))
    //{
    //  ShowErrorMessage("The specified SVM model directory does not exist.");
    //  return;
    //}
    ////--------------test subjects directory--------------------
    //if (testSubjectsDirectoryName->text().isEmpty())
    //{
    //  ShowErrorMessage("Please specify the directory of test subjects.");
    //  return;
    //}
    //if (!cbica::directoryExists(testSubjectsDirectoryName->text().toStdString()))
    //{
    //  ShowErrorMessage("The specified directory of test subjects does not exist.");
    //  return;
    //}
    ////--------------output directory--------------------
    //if (outputDirectoryName->text().isEmpty())
    //{
    //  ShowErrorMessage("Please specify the output directory.");
    //  return;
    //}
    //if (!cbica::directoryExists(outputDirectoryName->text().toStdString()))
    //{
    //  if (!cbica::createDirectory(outputDirectoryName->text().toStdString()))
    //  {
    //    ShowErrorMessage("Unable to create the output directory.");
    //    return;
    //  }
    //}
	//--------------function call--------------------
    emit EGFRvIIIPredictionOnExistingModel(svmModelFileName->text().toStdString(), testSubjectsDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
  else 
  {
    existingMaskDirectoryName->setEnabled(true);
    existingMasksButton->setEnabled(true);

	//--------------input subjects directory--------------------
    //if (existingMaskDirectoryName->text().isEmpty())
    //{
    //  ShowErrorMessage("Please specify the directory of training subjects.");
    //  return;
    //}
    //if (!cbica::directoryExists(existingMaskDirectoryName->text().toStdString()))
    //{
    //  ShowErrorMessage("The specified directory of training subjects does not exist.");
    //  return;
    //}
    ////--------------output directory--------------------
    //if (outputDirectoryName->text().isEmpty())
    //{
    //  ShowErrorMessage("Please specify the output directory.");
    //  return;
    //}
    //if (!cbica::directoryExists(outputDirectoryName->text().toStdString()))
    //{
    //  if (!cbica::createDirectory(outputDirectoryName->text().toStdString()))
    //  {
    //    ShowErrorMessage("Unable to create the output directory.");
    //    return;
    //  }
    //}
	//--------------function call--------------------
    emit PrepareNewEGFRvIIIPredictionModel(existingMaskDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
  this->close();
}

void fEGFRvIIIPredictor::OpenSVMModelFile()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    svmModelFileName->setText(directory);

  mInputPathName = directory;
}

void fEGFRvIIIPredictor::OpenExistingMasksDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    existingMaskDirectoryName->setText(directory);

  mInputPathName = directory;
}


void fEGFRvIIIPredictor::SelectOutputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);
}
void fEGFRvIIIPredictor::OpenTestSubjectsDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    testSubjectsDirectoryName->setText(directory);

  mInputPathName = directory;
}

void fEGFRvIIIPredictor::ExistingClassificationRadioButtonChecked()
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
void fEGFRvIIIPredictor::NewModelRadioButtonChecked()
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


void fEGFRvIIIPredictor::CheckForDisclaimer()
{
  QString volumeString;
  volumeString = "You are about to download a model trained on de novo glioblastoma cases.\n";
  volumeString += "Please note that this model was created following certain assumptions \n";
  volumeString += "(described in the paper below), can be used for research purposes only \n\n\n";
  volumeString += "H.Akbari, et al. In vivo evaluation of EGFRvIII mutation in primary\n";
  volumeString += "glioblastoma patients via complex multiparametric MRI signature, Neuro-Oncology. \n";
  volumeString += "20(8) : 1068–1079, 2018.";

  QMessageBox *box = new QMessageBox(QMessageBox::Question, "Disclaimer", volumeString, QMessageBox::Ok | QMessageBox::Cancel);
  //box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
  box->setWindowModality(Qt::NonModal);
  QCoreApplication::processEvents();

  if (box->exec() == QMessageBox::Ok)
  {
    cbica::Logging(loggerFile, m_trainedModelLink);

    //ShowErrorMessage("Starting download, may take a while, depending on your net bandwidth", this, "Downloading...");

    if /*(std::system((link).c_str()) != 0)*/ (!openLink(m_trainedModelLink))
    {
      ShowErrorMessage("CaPTk couldn't open the browser to download specified model.", this);
      return;
    }
    else
    {
      //std::string dataMessage = "Model has been saved to: " + captk_PretrainedFolder;
      //ShowMessage(dataMessage, this);
      return;
    }
  }
}