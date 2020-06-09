#include "fTrainingDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fTrainingSimulator::fTrainingSimulator()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(inputMaskButton, SIGNAL(clicked()), this, SLOT(OpenInputMaskImage()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(mSplitModelDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectSplitModelDirectory()));
  

}
fTrainingSimulator::~fTrainingSimulator()
{
}
void fTrainingSimulator::CancelButtonPressed()
{
  this->close();
}
void fTrainingSimulator::ConfirmButtonPressed()
{
  auto inputImageName_string = inputImageName->text().toStdString();
  auto outputImageName_string = outputImageName->text().toStdString();

  if ((inputImageName->text().isEmpty()))
  {
    ShowErrorMessage("Please select the features file.", this);
    return;
  }
  if ((inputMaskName->text().isEmpty()))
  {
    ShowErrorMessage("Please select the target file.", this);
    return;
  }
  if (!cbica::isFile(inputImageName->text().toStdString()))
  {
    ShowErrorMessage("The specified feature file does not exist.", this);
    return;
  }
  if (!cbica::isFile(inputMaskName->text().toStdString()))
  {
    ShowErrorMessage("The specified target file does not exist.", this);
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output directory.", this);
    return;
  }
  if (mLinearKernel->isChecked() == false && mRBFKernel->isChecked() == false)
  {
    ShowErrorMessage("Please select at least one of the given two options: Linear, RBF.");
    return;
  }
  if (mCrossValidation->isChecked() == false && mSplitTrainTest->isChecked() == false
      && mSplitTrain->isChecked() == false && mSplitTest->isChecked() == false)
  {
    ShowErrorMessage("Please select at least one of the given options: CrossValidation, TrainTest, Split Train, Split Test");
    return;
  }
  if (mCrossValidation->isChecked() == true && cvValue->text().isEmpty())
  {
    ShowErrorMessage("Please select the # of folds.");
    return;
  }
  if (mSplitTrainTest->isChecked() == true && ttValue->text().isEmpty())
  {
    ShowErrorMessage("Please select the size of training dataset.");
    return;
  }
  if (mSplitTest->isChecked() == true && mSplitModelDirectory->text().isEmpty())
  {
    ShowErrorMessage("Please provide a model directory to use for testing.");
    return;
  }


  int classifier = 0;
  int configuration = 0;
  int folds = 0;

  if (mLinearKernel->isChecked())
    classifier = 1;
  else
    classifier = 2;
  if (mCrossValidation->isChecked())
  {
    configuration = 1;
    folds = cvValue->text().toInt();
  }
  else if (mSplitTrainTest->isChecked())
  {
    configuration = 2;
    folds = ttValue->text().toInt();
  }
  else if (mSplitTrain->isChecked())
  {
    configuration = 3;
  }
  else
    configuration = 4;


  emit RunTrainingSimulation(mInputPathName.toStdString(), mInputMaskName.toStdString(), mOutputPathName.toStdString(),mModelDirectoryName.toStdString(), classifier, configuration, folds);
  this->close();
}



void fTrainingSimulator::OpenInputImage()
{
  QString extension_string = QString::fromStdString("CSV files: (*.csv)");
  auto inputImage = getExistingFile(this, mInputPathName,extension_string);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = inputImage;
}

void fTrainingSimulator::OpenInputMaskImage()
{
  QString extension_string = QString::fromStdString("CSV files: (*.csv)");
  auto inputImage = getExistingFile(this, mInputPathName, extension_string);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputMaskName->setText(inputImage);

  mInputMaskName = inputImage;
}


void fTrainingSimulator::SelectOutputImage()
{
  QString directory = getExistingDirectory(this, mOutputPathName);
  if (directory.isNull())
    return;
  else
    outputImageName->setText(directory);

  mOutputPathName = directory;
}

void fTrainingSimulator::SelectSplitModelDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    mSplitModelDirectory->setText(directory);

  mModelDirectoryName = directory;
}