#include "fRegistrationDialog.h"
#include "fProgressDialog.h"

fRegistrationDialog::fRegistrationDialog()
{
  setupUi(this);
  this->setWindowTitle("Register");
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));

  connect(fixedFileButton, SIGNAL(clicked()), this, SLOT(SelectFixedFile()));
  connect(movingFileButton1, SIGNAL(clicked()), this, SLOT(SelectMovingFile1()));
  connect(movingFileButton2, SIGNAL(clicked()), this, SLOT(SelectMovingFile2()));
  connect(movingFileButton3, SIGNAL(clicked()), this, SLOT(SelectMovingFile3()));
  connect(movingFileButton4, SIGNAL(clicked()), this, SLOT(SelectMovingFile4()));
  connect(movingFileButton5, SIGNAL(clicked()), this, SLOT(SelectMovingFile5()));

  connect(movingFileOutputButton1, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile1()));
  connect(movingFileOutputButton2, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile2()));
  connect(movingFileOutputButton3, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile3()));
  connect(movingFileOutputButton4, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile4()));
  connect(movingFileOutputButton5, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile5()));
}
fRegistrationDialog::~fRegistrationDialog()
{
}
//void fRegistrationDialog::ShowMessage(std::string msg)
//{
//	QMessageBox box(this);
//	box.setIcon(QMessageBox::Information);
//	box.addButton(QMessageBox::Ok);
//	box.setText(QString::fromStdString(msg));
//	box.setWindowTitle(tr("Missing Data"));
//	box.exec();
//	return;
//}

void fRegistrationDialog::SelectFixedFile()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  fixedFileName->setText(file);
}

void fRegistrationDialog::SelectMovingFile1()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  movingFileName1->setText(file);
}
void fRegistrationDialog::SelectMovingFile2()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  movingFileName2->setText(file);
}
void fRegistrationDialog::SelectMovingFile3()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  movingFileName3->setText(file);
}
void fRegistrationDialog::SelectMovingFile4()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  movingFileName4->setText(file);
}
void fRegistrationDialog::SelectMovingFile5()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  movingFileName5->setText(file);
}
void fRegistrationDialog::ConfirmButtonPressed()
{
  std::vector<std::string> inputfilenames;
  std::vector<std::string> outputfilenames;

  std::string mvFName1 = movingFileName1->text().toStdString(),
    mvFName2 = movingFileName2->text().toStdString(), mvFName3 = movingFileName3->text().toStdString(),
    mvFName4 = movingFileName4->text().toStdString(), mvFName5 = movingFileName5->text().toStdString();

  std::string mvOFName1 = movingFileOutputName1->text().toStdString(),
    mvOFName2 = movingFileOutputName2->text().toStdString(), mvOFName3 = movingFileOutputName3->text().toStdString(),
    mvOFName4 = movingFileOutputName4->text().toStdString(), mvOFName5 = movingFileOutputName5->text().toStdString();

  if (!mvFName1.empty())
    inputfilenames.push_back(mvFName1);
  if (!mvFName2.empty())
    inputfilenames.push_back(mvFName2);
  if (!mvFName3.empty())
    inputfilenames.push_back(mvFName3);
  if (!mvFName4.empty())
    inputfilenames.push_back(mvFName4);
  if (!mvFName5.empty())
    inputfilenames.push_back(mvFName5);

  if (!mvOFName1.empty())
    outputfilenames.push_back(mvOFName1);
  if (!mvOFName2.empty())
    outputfilenames.push_back(mvOFName2);
  if (!mvOFName3.empty())
    outputfilenames.push_back(mvOFName3);
  if (!mvOFName4.empty())
    outputfilenames.push_back(mvOFName4);
  if (!mvOFName5.empty())
    outputfilenames.push_back(mvOFName5);

  emit Registrationsignal(fixedFileName->text().toStdString(), inputfilenames, outputfilenames);
  this->close();
}
void fRegistrationDialog::CancelButtonPressed()
{
  this->close();
}

void fRegistrationDialog::SelectMovingOutputFile1()
{
  auto file = getSaveFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  movingFileOutputName1->setText(file);
}
void fRegistrationDialog::SelectMovingOutputFile2()
{
  auto file = getSaveFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  movingFileOutputName2->setText(file);
}
void fRegistrationDialog::SelectMovingOutputFile3()
{
  auto file = getSaveFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  movingFileOutputName3->setText(file);
}
void fRegistrationDialog::SelectMovingOutputFile4()
{
  auto file = getSaveFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  movingFileOutputName4->setText(file);
}
void fRegistrationDialog::SelectMovingOutputFile5()
{
  auto file = getSaveFile(this, mInputPathName);
  if (file.isEmpty())
  {
    return;
  }
  movingFileOutputName5->setText(file);
}