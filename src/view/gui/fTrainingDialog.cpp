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
  connect(inputFeaturesButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(inputTargetButton, SIGNAL(clicked()), this, SLOT(OpenInputMaskImage()));
  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(mSplitModelDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectSplitModelDirectory()));

  connect(mCrossValidation, SIGNAL(toggled(bool)), this, SLOT(CrossValidationRadioButtonChecked()));
  connect(mSplitTrain, SIGNAL(toggled(bool)), this, SLOT(TrainRadioButtonChecked()));
  connect(mSplitTest, SIGNAL(toggled(bool)), this, SLOT(TestRadioButtonChecked()));

  cvLabel->setEnabled(false);
  mSplitModelDirectoryLabel->setEnabled(false);
  cvValue->setEnabled(false);
  mSplitModelDirectory->setEnabled(false);
  mSplitModelDirectoryButton->setEnabled(false);
}
fTrainingSimulator::~fTrainingSimulator()
{
}
void fTrainingSimulator::CrossValidationRadioButtonChecked()
{
  if (mCrossValidation->isChecked())
  {
    cvLabel->setEnabled(true);
    cvValue->setEnabled(true);
    mSplitModelDirectoryLabel->setEnabled(false);
    mSplitModelDirectory->setEnabled(false);
    mSplitModelDirectoryButton->setEnabled(false);
  }
}
void fTrainingSimulator::SplitTrainRadioButtonChecked()
{
  if (mSplitTrain->isChecked())
  {
    cvLabel->setEnabled(false);
    cvValue->setEnabled(false);
    mSplitModelDirectoryLabel->setEnabled(false);
    mSplitModelDirectory->setEnabled(false);
    mSplitModelDirectoryButton->setEnabled(false);
  }
}
void fTrainingSimulator::SplitTestRadioButtonChecked()
{
  if (mSplitTest->isChecked())
  {
    cvLabel->setEnabled(false);
    cvValue->setEnabled(false);
    mSplitModelDirectoryLabel->setEnabled(true);
    mSplitModelDirectory->setEnabled(true);
    mSplitModelDirectoryButton->setEnabled(true);
  }
}
void fTrainingSimulator::CancelButtonPressed()
{
  this->close();
}
void fTrainingSimulator::ConfirmButtonPressed()
{
  if ((inputFeaturesName->text().isEmpty()))
  {
    ShowErrorMessage("Please select the features file.", this);
    return;
  }
  if (!cbica::isFile(inputFeaturesName->text().toStdString()))
  {
    ShowErrorMessage("The specified feature file does not exist.", this);
    return;
  }
  if (outputDirectoryName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output directory.", this);
    return;
  }

  if (mCrossValidation->isChecked() == false &&
      mSplitTrain->isChecked() == false &&
      mSplitTest->isChecked() == false)
  {
    ShowErrorMessage("Please select at least one of the given three options: CrossValidation, Train only, and Test only.");
    return;
  }

  if (mSplitTrain->isChecked() == true |
    mCrossValidation->isChecked() == true)
  {
    //error checks applied to the three configurations
    if (inputTargetName->text().isEmpty())
    {
      ShowErrorMessage("Please select the target file.", this);
      return;
    }
    if (!cbica::isFile(inputTargetName->text().toStdString()))
    {
      ShowErrorMessage("The specified target file does not exist.", this);
      return;
    }
    if (mLinearKernel->isChecked() == false && mRBFKernel->isChecked() == false)
    {
      ShowErrorMessage("Please select at least one of the given two classifiers: Linear, RBF.");
      return;
    }
    if (mSVMFFS->isChecked() == false && mESFS->isChecked() == false)
    {
      ShowErrorMessage("Please select at least one of the given two feature selection methods: SVM forward feature selection, Effect size feature selection.");
      return;
    }
    if (mResubstitution->isChecked() == false && mFiveFold->isChecked() == false)
    {
      ShowErrorMessage("Please select at least one of the given two cross-validation options: Resubstitution, 5-fold.");
      return;
    }

    //error checks applied to the individual configurations
    if (mCrossValidation->isChecked() == true && cvValue->text().isEmpty())
    {
      ShowErrorMessage("Please select the # of folds.");
      return;
    }
  }
  else if (mSplitTest->isChecked() == true)
  {
    if (mSplitModelDirectory->text().isEmpty())
    {
      ShowErrorMessage("The model directory field is empty.", this);
      return;
    }
    else if(!cbica::isDir(mSplitModelDirectory->text().toStdString()))
    {
      ShowErrorMessage("The specified model directory does not exist.", this);
      return;
    }
  }

  int classifierType = 1;
  int featureselectionType = 1;
  int optimizationType = 0;
  int crossvalidationType = 1;
  int foldType = 10;
  int confType = 1;

  std::string modelpath ="";

  if (mLinearKernel->isChecked())
    classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_LINEAR;
  else
    classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_RBF;

  if (mSVMFFS->isChecked())
    featureselectionType = CAPTK::FeatureSelectionType::FS_TYPE_FFS;
  else
    featureselectionType = CAPTK::FeatureSelectionType::FS_TYPE_ES;

  if (mResubstitution->isChecked())
    crossvalidationType = CAPTK::CrossValidationType::CV_TYPE_RESUBSTITUTION;
  else
    crossvalidationType = CAPTK::CrossValidationType::CV_TYPE_FiveFold;

  if (mOptimization->isChecked())
    optimizationType = CAPTK::OptimizationType::OPT_TYPE_ON;
  else
    optimizationType = CAPTK::OptimizationType::OPT_TYPE_OFF;

  if (mCrossValidation->isChecked())
  {
    confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_KFOLD_CV;
    foldType = cvValue->text().toInt();
  }
  else if (mSplitTrain->isChecked())
  {
    confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN;
  }
  else
  {
    confType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST;
    modelpath = mSplitModelDirectory->text().toStdString();
  }
  emit RunTrainingSimulation(mInputFeaturesName.toStdString(), mInputTargetName.toStdString(), mOutputPathName.toStdString(),mModelDirectoryName.toStdString(), classifierType, confType, foldType, featureselectionType, optimizationType,crossvalidationType);
  this->close();
}



void fTrainingSimulator::OpenInputImage()
{
  QString extension_string = QString::fromStdString("CSV files: (*.csv)");
  auto inputImage = getExistingFile(this, mInputFeaturesName,extension_string);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputFeaturesName->setText(inputImage);

  mInputFeaturesName = inputImage;
}

void fTrainingSimulator::OpenInputMaskImage()
{
  QString extension_string = QString::fromStdString("CSV files: (*.csv)");
  auto inputImage = getExistingFile(this, mInputTargetName, extension_string);
  if (inputImage.isNull() || inputImage.isEmpty())
    return;
  else
    inputTargetName->setText(inputImage);

  mInputTargetName = inputImage;
}


void fTrainingSimulator::SelectOutputImage()
{
  QString directory = getExistingDirectory(this, mOutputPathName);
  if (directory.isNull())
    return;
  else
    outputDirectoryName->setText(directory);

  mOutputPathName = directory;
}

void fTrainingSimulator::SelectSplitModelDirectory()
{
  QString directory = getExistingDirectory(this, mModelDirectoryName);
  if (directory.isNull())
    return;
  else
    mSplitModelDirectory->setText(directory);

  mModelDirectoryName = directory;
}