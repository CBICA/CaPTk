#include "fPseudoProgressionDialog.h"
#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "cbicaLogging.h"


fPseudoProgressionDialog::fPseudoProgressionDialog()
{
  this->setWindowTitle("Pseudo-Progression Estimator");
  //TBD read this from common place 
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  //this->setFixedWidth(400);
  //this->setFixedHeight(300);
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(existingMasksButton, SIGNAL(clicked()), this, SLOT(OpenExistingMasksDirectory()));
  connect(svmModelButton1, SIGNAL(clicked()), this, SLOT(OpenSVMModelFile1()));
  connect(svmModelButton2, SIGNAL(clicked()), this, SLOT(OpenSVMModelFile2()));
  connect(testSubjectsDirectoryButton, SIGNAL(clicked()), this, SLOT(OpenTestSubjectsDirectory()));
  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
  connect(rdExistingClassification, SIGNAL(toggled(bool)), this, SLOT(ExistingClassificationRadioButtonChecked()));
  connect(rdLoadedClassification, SIGNAL(toggled(bool)), this, SLOT(LoadedClassificationRadioButtonChecked()));
  connect(rdCreateModel, SIGNAL(toggled(bool)), this, SLOT(NewModelRadioButtonChecked()));
  connect(disclaimerButton, SIGNAL(clicked()), this, SLOT(CheckForDisclaimer()));

  //  connect(rdNewClassification, SIGNAL(toggled(bool)), this, SLOT(CurrentSubjectRadioButtonChecked()));

  //  rdNewClassification->setChecked(true);

  //rdNewClassification->hide();

  rdExistingClassification->setChecked(false);
  svmModelButton1->setEnabled(false);
  svmModelButton2->setEnabled(false);
  svmModelFileName1->setEnabled(false);
  svmModelFileName2->setEnabled(false);
  testSubjectsDirectoryButton->setEnabled(false);
  testSubjectsDirectoryName->setEnabled(false);
  existingMaskDirectoryName->setEnabled(false);
  existingMasksButton->setEnabled(false);
  cbT1Data->hide();
  cbDistanceData->hide();
  cbDTIData->hide();
  cbPerfData->hide();
  disclaimerLabel->setText("You can find a pretrained model, based on our study at PENN, ");

  rdLoadedClassification->hide();
  modelDirectoryLabel1->hide();
  //modelDirectoryLabel2->hide();
  svmModelFileName1->hide();
  svmModelButton1->hide();
}
fPseudoProgressionDialog::~fPseudoProgressionDialog()
{
}
void fPseudoProgressionDialog::CancelButtonPressed()
{
  this->close();
}
void fPseudoProgressionDialog::ConfirmButtonPressed()
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
  //------------------------loaded sbject checks------------------------------------------
  if (rdLoadedClassification->isChecked())
  {
    //if (svmModelFileName1->text().isEmpty())
    //{
    //  ShowErrorMessage("Please specify the directory of SVM model.");
    //  return;
    //}
    //if (cbica::isFile(svmModelFileName1->text().toStdString() + "/REC_SVM_Model.xml") && 
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/PSU_SVM_Model.xml") && 
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/PSU_ZScore_Mean.csv") && 
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/PSU_ZScore_Std.csv") &&
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/PSU_SelectedFeatures.csv") &&
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/REC_SelectedFeatures.csv") &&
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/Mean_Others.csv") && 
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/Mean_PERF.csv") && 
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/PCA_Others.csv") && 
    //  cbica::isFile(svmModelFileName1->text().toStdString() + "/PCA_PERF.csv"))
    //  emit SubjectBasedExistingPseudoprogressionEstimate(outputDirectoryName->text().toStdString(), svmModelFileName1->text().toStdString(), true, true, true, true);
    //else
    //{
    //  ShowErrorMessage("The specified directory does not have model files.");
    //  return;
    //}
    //this->close();
  }
  else if (rdExistingClassification->isChecked())
  {
    if (testSubjectsDirectoryName->text().isEmpty())
    {
      ShowErrorMessage("Please specify the directory of test subjects.");
      return;
    }
    if (svmModelFileName2->text().isEmpty())
    {
      ShowErrorMessage("Please specify the directory of SVM model.");
      return;
    }
    if (cbica::isFile(svmModelFileName2->text().toStdString() + "/REC_SVM_Model.xml") && 
      cbica::isFile(svmModelFileName2->text().toStdString() + "/PSU_SVM_Model.xml") && 
      cbica::isFile(svmModelFileName2->text().toStdString() + "/PSU_ZScore_Mean.csv") && 
      cbica::isFile(svmModelFileName2->text().toStdString() + "/PSU_ZScore_Std.csv") &&
      cbica::isFile(svmModelFileName2->text().toStdString() + "/PSU_SelectedFeatures.csv") &&
      cbica::isFile(svmModelFileName2->text().toStdString() + "/REC_SelectedFeatures.csv") &&
      cbica::isFile(svmModelFileName2->text().toStdString() + "/Mean_Others.csv") && 
      cbica::isFile(svmModelFileName2->text().toStdString() + "/Mean_PERF.csv") && 
      cbica::isFile(svmModelFileName2->text().toStdString() + "/PCA_Others.csv") && 
      cbica::isFile(svmModelFileName2->text().toStdString() + "/PCA_PERF.csv"))
      emit ExistingModelBasedPseudoprogressionEstimate(svmModelFileName2->text().toStdString(), testSubjectsDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString(), true, true, true, true);
    else
    {
      ShowErrorMessage("The specified directory does not have model files.");
      return;
    }
  }
  else
    emit TrainNewPseudoModel(existingMaskDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString(), true, true, true, true);
  this->close();
}

void fPseudoProgressionDialog::OpenSVMModelFile1()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    svmModelFileName1->setText(directory);

  mInputPathName = directory;
}
void fPseudoProgressionDialog::OpenSVMModelFile2()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    svmModelFileName2->setText(directory);

  mInputPathName = directory;
}



void fPseudoProgressionDialog::OpenExistingMasksDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    existingMaskDirectoryName->setText(directory);

  mInputPathName = directory;
}


void fPseudoProgressionDialog::SelectOutputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);
}
void fPseudoProgressionDialog::OpenTestSubjectsDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    testSubjectsDirectoryName->setText(directory);

  mInputPathName = directory;
}

void fPseudoProgressionDialog::ExistingClassificationRadioButtonChecked()
{
  if (rdExistingClassification->isChecked())
  {
    svmModelButton2->setEnabled(true);
    svmModelFileName2->setEnabled(true);
    svmModelButton1->setEnabled(false);
    svmModelFileName1->setEnabled(false);
    testSubjectsDirectoryButton->setEnabled(true);
    testSubjectsDirectoryName->setEnabled(true);
    existingMaskDirectoryName->setEnabled(false);
    existingMasksButton->setEnabled(false);
  }
}
void fPseudoProgressionDialog::LoadedClassificationRadioButtonChecked()
{
  if (rdLoadedClassification->isChecked())
  {
    svmModelButton1->setEnabled(true);
    svmModelFileName1->setEnabled(true);
    svmModelButton2->setEnabled(false);
    svmModelFileName2->setEnabled(false);
    testSubjectsDirectoryButton->setEnabled(false);
    testSubjectsDirectoryName->setEnabled(false);
    existingMaskDirectoryName->setEnabled(false);
    existingMasksButton->setEnabled(false);
  }
}


void fPseudoProgressionDialog::NewModelRadioButtonChecked()
{
  if (rdCreateModel->isChecked())
  {
    svmModelButton1->setEnabled(false);
    svmModelFileName1->setEnabled(false);
    svmModelButton2->setEnabled(false);
    svmModelFileName2->setEnabled(false);
    testSubjectsDirectoryButton->setEnabled(false);
    testSubjectsDirectoryName->setEnabled(false);
    existingMaskDirectoryName->setEnabled(true);
    existingMasksButton->setEnabled(true);
  }
}
void fPseudoProgressionDialog::CurrentSubjectRadioButtonChecked()
{
  //if (rdNewClassification->isChecked())
  //{
  //	svmModelButton->setEnabled(false);
  //	svmModelFileName->setEnabled(false);
  //	testSubjectsDirectoryButton->setEnabled(false);
  //	testSubjectsDirectoryName->setEnabled(false);
  //	existingMaskDirectoryName->setEnabled(false);
  //	existingMasksButton->setEnabled(false);
  //}
}

void fPseudoProgressionDialog::CheckForDisclaimer()
{
  QString volumeString;
  volumeString = "You are about to download a model trained on de novo glioblastoma cases.\n";
  volumeString += "Please note that this model was created following certain assumptions \n";
  volumeString += "(described in the paper below), can be used for research purposes only \n\n\n";
  volumeString += "H.Akbari, et al. Quantitative radiomics and machine learning to distinguish\n";
  volumeString += "true progression from pseudoprogression in patients with GBM, \n";
  volumeString += "ASNR 56th Annual Meeting, 2018.\n";

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
      //ShowMessage(dataMessage, this, "Saved");
      return;
    }
  }
}
