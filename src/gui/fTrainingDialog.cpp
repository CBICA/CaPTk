#include "fTrainingDialog.h"
#include "fProgressDialog.h"

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

  if (inputFoldsName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the number of folds.", this);
    return;
  }
  if (mLinearKernel->isChecked() == false && mRBFKernel->isChecked() == false)
  {
    ShowErrorMessage("Please select at least one of the given two options: Linear, RBF.");
    return;
  }
  if(mLinearKernel->isChecked())
    emit RunTrainingSimulation(mInputPathName.toStdString(), mInputMaskName.toStdString(), mOutputPathName.toStdString(),1, inputFoldsName->text().toInt());

  if (mRBFKernel->isChecked())
    emit RunTrainingSimulation(mInputPathName.toStdString(), mInputMaskName.toStdString(), mOutputPathName.toStdString(), 2, inputFoldsName->text().toInt());


  this->close();
}



void fTrainingSimulator::OpenInputImage()
{
  QString extension_string = QString::fromStdString("CSV files: (*.csv)");
  auto inputImage = getExistingFile(this, mInputPathName,extension_string);
  if (inputImage.isNull())
    return;
  else
    inputImageName->setText(inputImage);

  mInputPathName = inputImage;
}

void fTrainingSimulator::OpenInputMaskImage()
{
  QString extension_string = QString::fromStdString("CSV files: (*.csv)");
  auto inputImage = getExistingFile(this, mInputPathName, extension_string);
  if (inputImage.isNull())
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
