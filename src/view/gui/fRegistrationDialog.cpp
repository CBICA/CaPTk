#include "fRegistrationDialog.h"
#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "itkFileTools.h"
#include <sys/types.h>
#include <sys/stat.h>
#ifndef WIN32
#include <unistd.h>
#endif

#ifdef WIN32
#define stat _stat
#endif

fRegistrationDialog::fRegistrationDialog()
{
  setupUi(this);
  this->setWindowTitle("Register");
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(resetButton, SIGNAL(clicked()), this, SLOT(ResetButtonPressed()));

  connect(fixedFileButton, SIGNAL(clicked()), this, SLOT(SelectFixedFile()));
  connect(movingFileButton1, SIGNAL(clicked()), this, SLOT(SelectMovingFile1()));
  connect(movingFileButton2, SIGNAL(clicked()), this, SLOT(SelectMovingFile2()));
  connect(movingFileButton3, SIGNAL(clicked()), this, SLOT(SelectMovingFile3()));
  connect(movingFileButton4, SIGNAL(clicked()), this, SLOT(SelectMovingFile4()));
  connect(movingFileButton5, SIGNAL(clicked()), this, SLOT(SelectMovingFile5()));


  options_AFFINE_selected->setChecked(true);
  //nccRadii->setDisabled(true);
  confirmButton->setText(QString("Register"));

  movingFileLabel2->setHidden(true);
  movingFileName2->setHidden(true);
  movingFileButton2->setHidden(true);
  movingFileOutputLabel2->setHidden(true);
  movingFileOutputName2->setHidden(true);
  movingFileOutputButton2->setHidden(true);

  movingFileLabel3->setHidden(true);
  movingFileName3->setHidden(true);
  movingFileButton3->setHidden(true);
  movingFileOutputLabel3->setHidden(true);
  movingFileOutputName3->setHidden(true);
  movingFileOutputButton3->setHidden(true);

  movingFileLabel4->setHidden(true);
  movingFileName4->setHidden(true);
  movingFileButton4->setHidden(true);
  movingFileOutputLabel4->setHidden(true);
  movingFileOutputName4->setHidden(true);
  movingFileOutputButton4->setHidden(true);

  movingFileLabel5->setHidden(true);
  movingFileName5->setHidden(true);
  movingFileButton5->setHidden(true);
  movingFileOutputLabel5->setHidden(true);
  movingFileOutputName5->setHidden(true);
  movingFileOutputButton5->setHidden(true);

  matrixLabel2->setHidden(true);
  matrixName2->setHidden(true);
  matrixButton2->setHidden(true);

  matrixLabel3->setHidden(true);
  matrixName3->setHidden(true);
  matrixButton3->setHidden(true);

  matrixLabel4->setHidden(true);
  matrixName4->setHidden(true);
  matrixButton4->setHidden(true);

  matrixLabel5->setHidden(true);
  matrixName5->setHidden(true);
  matrixButton5->setHidden(true);

  matrixRadioButton->setChecked(true);

  connect(matrixRadioButton, SIGNAL(clicked(bool)), this, SLOT(SelectGenerateMatrix(bool)));
  connect(matrixButton1, SIGNAL(clicked()), this, SLOT(SelectMatrixFile1()));
  connect(matrixButton2, SIGNAL(clicked()), this, SLOT(SelectMatrixFile2()));
  connect(matrixButton3, SIGNAL(clicked()), this, SLOT(SelectMatrixFile3()));
  connect(matrixButton4, SIGNAL(clicked()), this, SLOT(SelectMatrixFile4()));
  connect(matrixButton5, SIGNAL(clicked()), this, SLOT(SelectMatrixFile5()));
  connect(movingFileOutputButton1, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile1()));
  connect(movingFileOutputButton2, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile2()));
  connect(movingFileOutputButton3, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile3()));
  connect(movingFileOutputButton4, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile4()));
  connect(movingFileOutputButton5, SIGNAL(clicked()), this, SLOT(SelectMovingOutputFile5()));

  connect(addMoreButton, SIGNAL(clicked()), this, SLOT(addMoreImages()));

  connect(options_AFFINE_selected, SIGNAL(toggled(bool)), this, SLOT(SelectedAffineMode()));
  connect(options_RIGID_selected, SIGNAL(toggled(bool)), this, SLOT(SelectedRigidMode()));
  connect(options_DEFORMABLE_selected, SIGNAL(toggled(bool)), this, SLOT(SelectedDeformMode()));

  connect(options_MetricSelector, SIGNAL(currentIndexChanged(int)), this, SLOT(SelectedMetric(int)));
  connect(nccRadii, SIGNAL(textChanged(QString)), this, SLOT(setRadii(QString)));
  connect(iterations, SIGNAL(textChanged(QString)), this, SLOT(getIterations(QString)));

}

fRegistrationDialog::~fRegistrationDialog()
{
  /*
  QMessageBox box(this);
  box.setIcon(QMessageBox::Information);
  */

}

/*
void fRegistrationDialog::ShowMessage(std::string msg)
{
  QMessageBox box(this);
  box.setIcon(QMessageBox::Information);
  box.addButton(QMessageBox::Ok);
  box.setText(QString::fromStdString(msg));
  box.setWindowTitle(tr("Missing Data"));
  box.exec();
  return;
}
*/

void fRegistrationDialog::SelectFixedFile()
{
  auto file = getExistingFile(this, mInputPathName);
  if (file.isEmpty() || file.isNull())
  {
    return;
  }
  mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
  fixedFileName->setText(file);
}

void fRegistrationDialog::populateFields(std::string &configfilePath, std::string &matfilePath) {

  if (cbica::fileExists(configfilePath) && cbica::fileExists(matfilePath)) {

    //ShowErrorMessage("A transformation matrix was already present for this moving image. Parameters have been pre-populated for reference.");
    //auto filename = configfilePath;
    //struct stat result;
    //if (stat(filename.c_str(), &result) == 0)
    //{
    //    auto mod_time = result.st_mtime;
    //    char timebuf[26];
    //    errno_t err;
    //    // Get UNIX-style time and display as number and string.
    //    time(&mod_time);
    //    ("Time in seconds since UTC 1/1/70:\t%lld\n", (long long)mod_time);
    //    err = ctime_s(timebuf, 26, &mod_time);
    //    if (err)
    //    {
    //        printf("ctime_s failed due to an invalid argument.");
    //        exit(1);
    //    }
    //    std::string output = ("UNIX time and date:\t\t\t%s", timebuf);
    //    //ShowErrorMessage(output);
    //}

    std::string timestamp;
    std::vector<std::string> paramValues;
    std::ifstream inputFile(configfilePath.c_str());

    if (!inputFile)
    {
      std::cerr << "File '" << configfilePath << "' not found.\n";
      exit(EXIT_FAILURE);
    }
    std::string line;
    const char comma = ',';

    while (std::getline(inputFile, line, comma))
    {
      paramValues.push_back(line);
    }

    inputFile.close();

    /*if (paramValues.size() == 6)
        timestamp = paramValues[5];

    else if (paramValues.size() == 5)
        timestamp = paramValues[4];
    else
        timestamp = "No timestamp found!!";

    long int i = atol(timestamp.c_str());
    time_t newTime = static_cast<time_t> (i);
    auto mainTime = ctime(&newTime);
    ShowErrorMessage(mainTime);*/

    QMessageBox::StandardButton reply;

    reply = QMessageBox::question(this, "Warning!!", "A transformation matrix is already present from for this moving image. Do you want to re-use those parameters?",
      QMessageBox::Yes | QMessageBox::No);

    if (reply == QMessageBox::Yes) {
      int index = 0;
      if (paramValues[1] == "NMI") {
        options_MetricSelector->setCurrentIndex(0);
      }
      else if (paramValues[1] == "MI") {
        options_MetricSelector->setCurrentIndex(1);
      }
      else if (paramValues[1] == "NCC") {
        options_MetricSelector->setCurrentIndex(2);
        //Patch radii to specify which area to focus registration on
        nccRadii->setText(QString::fromStdString(paramValues[2]));
      }
      else if (paramValues[1] == "SSD") {
        options_MetricSelector->setCurrentIndex(3);
      }

      fixedFileName->setText(QString::fromStdString(paramValues[0]));

      if (paramValues.size() == 6) {
        if (paramValues[3] == "Affine") {
          options_AFFINE_selected->setChecked(true);
        }
        else {
          options_RIGID_selected->setChecked(true);
        }
        //Number of high, medium and low resolution iterations
        iterations->setText(QString::fromStdString(paramValues[4]));
      }
      else {
        if (paramValues[2] == "Affine") {
          options_AFFINE_selected->setChecked(true);
        }
        else {
          options_RIGID_selected->setChecked(true);
        }
        iterations->setText(QString::fromStdString(paramValues[3]));
      }
    }
  }
  else {
    return;
  }
}

void fRegistrationDialog::SelectMovingFile1()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") 
  {
    auto file = getExistingFile(this, mInputPathName);
    if (file.isEmpty() || file.isNull())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    QString fileName = itksys::SystemTools::GetFilenameWithoutExtension(file.toStdString()).c_str();
    QString extn = itksys::SystemTools::GetFilenameExtension(file.toStdString()).c_str();
    QString fixedfileName = itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str();

    movingFileName1->setText(file);
    if (!matrixRadioButton->isChecked())
      matrixName1->setText(mInputPathName + "/" + fileName + "_reg_to_remove_" + fixedfileName + matExtn);
    else
      matrixName1->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + matExtn);

    movingFileOutputName1->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + extn);

    std::string configFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".txt";

    std::string matFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".mat";
    populateFields(configFile, matFile);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMovingFile2()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") 
  {
    auto file = getExistingFile(this, mInputPathName);
    if (file.isEmpty() || file.isNull())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    QString fileName = itksys::SystemTools::GetFilenameWithoutExtension(file.toStdString()).c_str();
    QString extn = itksys::SystemTools::GetFilenameExtension(file.toStdString()).c_str();
    QString fixedfileName = itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str();

    movingFileName2->setText(file);

    if (!matrixRadioButton->isChecked())
      matrixName2->setText(mInputPathName + "/" + fileName + "_reg_to_remove_" + fixedfileName + matExtn);
    else
      matrixName2->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + matExtn);

    movingFileOutputName2->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + extn);

    std::string configFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".txt";

    std::string matFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".mat";

    populateFields(configFile, matFile);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMovingFile3()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") 
  {
    auto file = getExistingFile(this, mInputPathName);
    if (file.isEmpty() || file.isNull())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    QString fileName = itksys::SystemTools::GetFilenameWithoutExtension(file.toStdString()).c_str();
    QString extn = itksys::SystemTools::GetFilenameExtension(file.toStdString()).c_str();
    QString fixedfileName = itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str();

    movingFileName3->setText(file);

    if (!matrixRadioButton->isChecked())
      matrixName3->setText(mInputPathName + "/" + fileName + "_reg_to_remove_" + fixedfileName + matExtn);
    else
      matrixName3->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + matExtn);

    movingFileOutputName3->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + extn);

    std::string configFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".txt";

    std::string matFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".mat";
    populateFields(configFile, matFile);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMovingFile4()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "")
  {
    auto file = getExistingFile(this, mInputPathName);
    if (file.isEmpty() || file.isNull())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    QString fileName = itksys::SystemTools::GetFilenameWithoutExtension(file.toStdString()).c_str();
    QString extn = itksys::SystemTools::GetFilenameExtension(file.toStdString()).c_str();
    QString fixedfileName = itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str();

    movingFileName4->setText(file);

    if (!matrixRadioButton->isChecked())
      matrixName4->setText(mInputPathName + "/" + fileName + "_reg_to_remove_" + fixedfileName + matExtn);
    else
      matrixName4->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + matExtn);

    movingFileOutputName4->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + extn);

    std::string configFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".txt";

    std::string matFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".mat";
    populateFields(configFile, matFile);
  }
  else
  {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMovingFile5()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") 
  {
    auto file = getExistingFile(this, mInputPathName);
    if (file.isEmpty() || file.isNull())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    QString fileName = itksys::SystemTools::GetFilenameWithoutExtension(file.toStdString()).c_str();
    QString extn = itksys::SystemTools::GetFilenameExtension(file.toStdString()).c_str();
    QString fixedfileName = itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str();

    movingFileName5->setText(file);

    if (!matrixRadioButton->isChecked())
      matrixName5->setText(mInputPathName + "/" + fileName + "_reg_to_remove_" + fixedfileName + matExtn);
    else
      matrixName5->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + matExtn);

    movingFileOutputName5->setText(mInputPathName + "/" + fileName + "_reg_to_" + fixedfileName + extn);

    std::string configFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".txt";

    std::string matFile = mInputPathName.toStdString() + "/" + fileName.toStdString() + "_reg_to_" + itksys::SystemTools::GetFilenameWithoutExtension(fixed).c_str() + ".mat";
    populateFields(configFile, matFile);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMatrixFile1()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName, "", "Matrix (*.mat)");
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    matrixName1->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMatrixFile2()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName, "", "Matrix (*.mat)");
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    matrixName2->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMatrixFile3()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName, "", "Matrix (*.mat)");
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    matrixName3->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMatrixFile4()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName, "", "Matrix (*.mat)");
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    matrixName4->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}


void fRegistrationDialog::SelectMatrixFile5()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName, "", "Matrix (*.mat)");
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    matrixName5->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMovingOutputFile1()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName);
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    movingFileOutputName1->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}
void fRegistrationDialog::SelectMovingOutputFile2()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName);
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    movingFileOutputName2->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}
void fRegistrationDialog::SelectMovingOutputFile3()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName);
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    movingFileOutputName3->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}
void fRegistrationDialog::SelectMovingOutputFile4()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName);
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    movingFileOutputName4->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectMovingOutputFile5()
{
  auto fixed = fixedFileName->text().toStdString();
  if (fixed != "") {
    auto file = getSaveFile(this, mInputPathName);
    if (file.isEmpty())
    {
      return;
    }
    mInputPathName = itksys::SystemTools::GetFilenamePath(file.toStdString()).c_str();
    movingFileOutputName5->setText(file);
  }
  else {
    ShowErrorMessage("Select fixed image first!", this, "Warning!");;
    fixedFileName->setFocus();
  }
}

void fRegistrationDialog::SelectedAffineMode()
{
  options_RIGID_selected->setChecked(false);
  options_DEFORMABLE_selected->setChecked(false);
  affineMode = true;
  rigidMode = false;
  deformMode = false;
}

void fRegistrationDialog::SelectedRigidMode()
{
  options_AFFINE_selected->setChecked(false);
  options_DEFORMABLE_selected->setChecked(false);
  affineMode = false;
  rigidMode = true;
  deformMode = false;
}

void fRegistrationDialog::SelectedDeformMode()
{
  options_AFFINE_selected->setChecked(false);
  options_RIGID_selected->setChecked(false);
  affineMode = false;
  rigidMode = false;
  deformMode = true;
}

void fRegistrationDialog::SelectedMetric(int index)
{
  switch (options_MetricSelector->currentIndex())
  {
  case 0:
  {
    metric = "NMI";
    options_NCC_radii->setDisabled(true);
    nccRadii->setDisabled(true);
    break;
  }
  case 1:
  {
    metric = "MI";
    options_NCC_radii->setDisabled(true);
    nccRadii->setDisabled(true);
    break;
  }
  case 2:
  {
    metric = "NCC";
    radii = true;
    options_NCC_radii->setDisabled(false);
    nccRadii->setDisabled(false);
    break;
  }
  case 3:
  {
    metric = "SSD";
    options_NCC_radii->setDisabled(true);
    nccRadii->setDisabled(true);
    break;
  }
  default:
  {
    metric = "NMI";
    options_NCC_radii->setDisabled(true);
    nccRadii->setDisabled(true);
    break;
  }
  }
}

void fRegistrationDialog::addMoreImages()
{
  clicked++;

  if (clicked == 1) {
    movingFileLabel2->setHidden(false);
    movingFileName2->setHidden(false);
    movingFileButton2->setHidden(false);
    matrixLabel2->setHidden(false);
    matrixName2->setHidden(false);
    matrixButton2->setHidden(false);
    movingFileOutputLabel2->setHidden(false);
    movingFileOutputName2->setHidden(false);
    movingFileOutputButton2->setHidden(false);
  }
  else if (clicked == 2) {
    movingFileLabel3->setHidden(false);
    movingFileName3->setHidden(false);
    movingFileButton3->setHidden(false);
    matrixLabel3->setHidden(false);
    matrixName3->setHidden(false);
    matrixButton3->setHidden(false);
    movingFileOutputLabel3->setHidden(false);
    movingFileOutputName3->setHidden(false);
    movingFileOutputButton3->setHidden(false);
  }
  else if (clicked == 3) {
    movingFileLabel4->setHidden(false);
    movingFileName4->setHidden(false);
    movingFileButton4->setHidden(false);
    matrixLabel4->setHidden(false);
    matrixName4->setHidden(false);
    matrixButton4->setHidden(false);
    movingFileOutputLabel4->setHidden(false);
    movingFileOutputName4->setHidden(false);
    movingFileOutputButton4->setHidden(false);
  }
  else if (clicked == 4) {
    movingFileLabel5->setHidden(false);
    movingFileName5->setHidden(false);
    movingFileButton5->setHidden(false);
    matrixLabel5->setHidden(false);
    matrixName5->setHidden(false);
    matrixButton5->setHidden(false);
    movingFileOutputLabel5->setHidden(false);
    movingFileOutputName5->setHidden(false);
    movingFileOutputButton5->setHidden(false);
  }
}

void fRegistrationDialog::SelectGenerateMatrix(bool checked)
{
  if (checked)
    matrixGroupBox->setHidden(false);
  else
    matrixGroupBox->setHidden(true);
  //QLayout * layout = this->layout();
  //layout->setSizeConstraint(QLayout::SetFixedSize);
  clicked = 0;
}


void fRegistrationDialog::ResetButtonPressed() {
  //QLayout * layout = this->layout();
  //layout->setSizeConstraint(QLayout::SetFixedSize);
  clicked = 0;
  options_AFFINE_selected->setChecked(true);
  options_MetricSelector->setCurrentIndex(0);
  nccRadii->setText("5x5x5");
  iterations->setText("100x50x5");

  movingFileName1->setText("");
  movingFileName2->setText("");
  movingFileName3->setText("");
  movingFileName4->setText("");
  movingFileName5->setText("");
  movingFileOutputName1->setText("");
  movingFileOutputName2->setText("");
  movingFileOutputName3->setText("");
  movingFileOutputName4->setText("");
  movingFileOutputName5->setText("");
  fixedFileName->setText("");
  matrixName1->setText("");
  matrixName2->setText("");
  matrixName3->setText("");
  matrixName4->setText("");
  matrixName5->setText("");


  movingFileLabel2->setHidden(true);
  movingFileName2->setHidden(true);
  movingFileButton2->setHidden(true);
  movingFileOutputLabel2->setHidden(true);
  movingFileOutputName2->setHidden(true);
  movingFileOutputButton2->setHidden(true);
  matrixLabel2->setHidden(true);
  matrixName2->setHidden(true);
  matrixButton2->setHidden(true);

  movingFileLabel3->setHidden(true);
  movingFileName3->setHidden(true);
  movingFileButton3->setHidden(true);
  movingFileOutputLabel3->setHidden(true);
  movingFileOutputName3->setHidden(true);
  movingFileOutputButton3->setHidden(true);
  matrixLabel3->setHidden(true);
  matrixName3->setHidden(true);
  matrixButton3->setHidden(true);

  movingFileLabel4->setHidden(true);
  movingFileName4->setHidden(true);
  movingFileButton4->setHidden(true);
  movingFileOutputLabel4->setHidden(true);
  movingFileOutputName4->setHidden(true);
  movingFileOutputButton4->setHidden(true);
  matrixLabel4->setHidden(true);
  matrixName4->setHidden(true);
  matrixButton4->setHidden(true);

  movingFileLabel5->setHidden(true);
  movingFileName5->setHidden(true);
  movingFileButton5->setHidden(true);
  movingFileOutputLabel5->setHidden(true);
  movingFileOutputName5->setHidden(true);
  movingFileOutputButton5->setHidden(true);
  matrixLabel5->setHidden(true);
  matrixName5->setHidden(true);
  matrixButton5->setHidden(true);
}


void fRegistrationDialog::setRadii(const QString nccRadii)
{
  radius = nccRadii.toStdString();
}

void fRegistrationDialog::getIterations(const QString iterations)
{
  m_iterations = iterations.toStdString();
}

void fRegistrationDialog::ConfirmButtonPressed()
{
  std::vector<std::string> inputfilenames;
  std::vector<std::string> outputfilenames;
  std::vector<std::string> matrixfilenames;

  std::string affineMatrix;
  std::string warpImage;
  std::string weights;

  std::string mvFName1 = movingFileName1->text().toStdString(),
    mvFName2 = movingFileName2->text().toStdString(), mvFName3 = movingFileName3->text().toStdString(),
    mvFName4 = movingFileName4->text().toStdString(), mvFName5 = movingFileName5->text().toStdString();

  std::string mvOFName1 = movingFileOutputName1->text().toStdString(),
    mvOFName2 = movingFileOutputName2->text().toStdString(), mvOFName3 = movingFileOutputName3->text().toStdString(),
    mvOFName4 = movingFileOutputName4->text().toStdString(), mvOFName5 = movingFileOutputName5->text().toStdString();

  std::string matName1 = matrixName1->text().toStdString(),
    matName2 = matrixName2->text().toStdString(), matName3 = matrixName3->text().toStdString(),
    matName4 = matrixName4->text().toStdString(), matName5 = matrixName5->text().toStdString();

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
  if (!matName1.empty())
    matrixfilenames.push_back(matName1);
  if (!matName2.empty())
    matrixfilenames.push_back(matName2);
  if (!matName3.empty())
    matrixfilenames.push_back(matName3);
  if (!matName4.empty())
    matrixfilenames.push_back(matName4);
  if (!matName5.empty())
    matrixfilenames.push_back(matName5);

  if ((fixedFileName->text().toStdString()).empty()) {
    ShowErrorMessage("Fixed file cannot be empty");
    fixedFileName->setFocus();
  }

  if ((fixedFileName->text().toStdString()).empty() || inputfilenames.size() != outputfilenames.size() || inputfilenames.size() != matrixfilenames.size()
    || outputfilenames.size() != matrixfilenames.size() || inputfilenames.size() < 1 || matrixfilenames.size() < 1 || outputfilenames.size() < 1) {
    ShowErrorMessage("Number of input, output and matrix file names do not match. Please check the fields.");
  }

  else 
  {
    emit RegistrationSignal(fixedFileName->text().toStdString(), inputfilenames, 
      outputfilenames, matrixfilenames, 
      metric, rigidMode, affineMode, deformMode, radius, m_iterations);
    this->close();
  }
}

void fRegistrationDialog::CancelButtonPressed()
{
  this->close();
}

//void fRegistrationDialog::SelectedRegistration()
//{
//    confirmButton->setObjectName(QString::fromUtf8("Generate"));
//    option_transformationMode->setChecked(false);
//    registrationMode = true;
//    registrationGroupBox->setDisabled(false);
//    addMoreButton->setHidden(false);
//    //movingFileOutputLabel1->setText("Select Matrix:");
//    confirmButton->setText(QString("Register"));
//    outputGroupBox->setHidden(true);
//    modeGroupBox->setDisabled(false);
//    registrationGroupBox->setDisabled(false);
//    /* options_MetricSelector->setHidden(false);
//    options_MetricSelector_label->setHidden(false);
//    options_AFFINE_selected->setHidden(false);
//    options_RIGID_selected->setHidden(false);
//    options_registration_label->setHidden(false);*/
//
//}
//
//void fRegistrationDialog::SelectedTransformation()
//{
//    confirmButton->setObjectName(QString::fromUtf8("Apply"));
//    option_registrationMode->setChecked(false);
//    registrationMode = false;
//    registrationGroupBox->setDisabled(true);
//    addMoreButton->setHidden(false);
//    confirmButton->setText(QString("Apply"));
//    //movingFileOutputLabel1->setText("Output Image 1:");
//    outputGroupBox->setHidden(false);
//    modeGroupBox->setDisabled(true);
//    registrationGroupBox->setDisabled(true);
//    /*options_MetricSelector->setHidden(true);
//    options_MetricSelector_label->setHidden(true);
//    options_registration_label->setHidden(true);
//    options_AFFINE_selected->setHidden(true);
//    options_RIGID_selected->setHidden(true);    */
//
//}
