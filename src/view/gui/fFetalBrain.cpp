#include "fFetalBrain.h"
#include "FetalBrain.h"
#include "CaPTkGUIUtils.h"

fFetalBrain::fFetalBrain()
{
  this->setWindowTitle("Fetal Ventricomegaly");
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  ////TBD read this from common place 
  setupUi(this);
  
  connect(Segment, SIGNAL(clicked()), this, SLOT(SkullSegmentation()));
  connect(Predict, SIGNAL(clicked()), this, SLOT(Prediction()));

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));

  connect(rdExistingClassification, SIGNAL(toggled(bool)), this, SLOT(SingleSubject()));
  connect(rdCreateModel, SIGNAL(toggled(bool)), this, SLOT(TrainingRadioButtonChecked()));

  connect(TrainingSubDirButton, SIGNAL(clicked()), this, SLOT(OpenTrainSubjectDirectory()));

  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputDirectory()));
  std::string m_tempFolderLocation=  loggerFolder + "tmp_" + cbica::getCurrentProcessID();
  outputDirectoryName->setText(QString::fromStdString(m_tempFolderLocation));
}

fFetalBrain::~fFetalBrain()
{
}

void fFetalBrain::SkullSegmentation()
{
  if (rdExistingClassification->isChecked())
  {
    emit skullstripfun();
    Predict->setEnabled(true);
  }
}

void fFetalBrain::Prediction()
{

  emit drawlinear();
 
}


void fFetalBrain::CancelButtonPressed()
{
  this->close();
}
void fFetalBrain::ConfirmButtonPressed()
{
  if (rdCreateModel->isChecked())
  {
    emit TrainNewFetalModel(trainDirectoryName->text().toStdString(), outputDirectoryName->text().toStdString());
  }
}

void fFetalBrain::OpenSVMModelFile()
{
 
QString filename = QFileDialog::getOpenFileName(this,
    tr("Open Model File"), QString(mInputPathName.c_str()), tr("XML files (*.xml)"));
QFileInfo check_file(filename);
if (check_file.exists() && check_file.isFile())
{
 // svmModelFileName->setText(filename);
}
else{
  ShowErrorMessage("File doesn't exist");
  return;
}
  
}


void fFetalBrain::SelectOutputDirectory()
{
  std::string root_directory;
  //int imagetype;
  QString directory = QFileDialog::getExistingDirectory(this, tr("Open Directory"), QString(mInputPathName.c_str()), QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);
}
void fFetalBrain::OpenTestSubjectsDirectory()
{
  std::string root_directory;
  //int imagetype;
  QString directory = QFileDialog::getExistingDirectory(this, tr("Open Directory"), QString(mInputPathName.c_str()), QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  if (directory.isNull())
    return;
  else
  { }
  //  testSubjectsDirectoryName->setText(directory);
}

std::string fFetalBrain::OpenTrainSubjectDirectory()
{
  QString directory = QFileDialog::getExistingDirectory(this, tr("Open Directory"), QString(mInputPathName.c_str()), QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  if (directory.isNull())
    return "";
  else
    trainDirectoryName->setText(directory);
  return trainDirectoryName->accessibleName().toUtf8().constData();
}

void fFetalBrain::SingleSubject()
{
  if (rdExistingClassification->isChecked())
  {
    Segment->setEnabled(true);
    Predict->setEnabled(false);
    trainDirectoryLabel->setEnabled(false);
    trainDirectoryName->setEnabled(false);
    TrainingSubDirButton->setEnabled(false);
    outputDirectoryButton->setEnabled(false);
    outputDirectoryLabel->setEnabled(false);
    outputDirectoryName->setEnabled(false);

    rdCreateModel->setChecked(false);
    confirmButton->setEnabled(false);
  }
}



void fFetalBrain::TrainingRadioButtonChecked()
{
  if (rdCreateModel->isChecked())
  {
    trainDirectoryLabel->setEnabled(true);
    trainDirectoryName->setEnabled(true);
    TrainingSubDirButton->setEnabled(true);
    outputDirectoryButton->setEnabled(true);
    outputDirectoryLabel->setEnabled(true);
    outputDirectoryName->setEnabled(true);
    confirmButton->setEnabled(true);

    rdExistingClassification->setChecked(false);

    Segment->setEnabled(false);
    Predict->setEnabled(false);
  }
}
