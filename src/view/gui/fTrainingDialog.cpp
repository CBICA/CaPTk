#include "fTrainingDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"
#include "TrainingModule.h"

fTrainingSimulator::fTrainingSimulator()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  // This can be re-enabled once Qt is fully enabled for TrainingModule
  //connect(&m_trainingSimulator, SIGNAL(updateProgress(int)), this, SLOT(onProgressUpdate(int)));
  //connect(&m_trainingSimulator, SIGNAL(done(TrainingModuleResult)), this, SLOT(onJobDone(TrainingModuleResult)));

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputFeaturesButton, SIGNAL(clicked()), this, SLOT(OpenInputImage()));
  connect(inputTargetButton, SIGNAL(clicked()), this, SLOT(OpenInputMaskImage()));
  connect(outputDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(mSplitModelDirectoryButton, SIGNAL(clicked()), this, SLOT(SelectSplitModelDirectory()));

  connect(mCrossValidation, SIGNAL(toggled(bool)), this, SLOT(CrossValidationRadioButtonChecked()));
  connect(mSplitTrain, SIGNAL(toggled(bool)), this, SLOT(SplitTrainRadioButtonChecked()));
  connect(mSplitTest, SIGNAL(toggled(bool)), this, SLOT(SplitTestRadioButtonChecked()));
  connect(mOptimization, SIGNAL(toggled(bool)), this, SLOT(OptimizationToggled(bool)));

  connect(mLinearKernel, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));
  connect(mRBFKernel, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));
  connect(mPolynomialKernel, SIGNAL(clicked(bool)), this, SLOT(SVMToggled(bool)));
  connect(mSigmoidKernel, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));
  connect(mChiSquaredKernel, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));
  connect(mIntersectionKernel, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));
  
  // Just included to help clear options from the GUI, this needs to be reorganized
  connect(mSGDSVM, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));
  connect(mBoostedTrees, SIGNAL(toggled(bool)), this, SLOT(SVMToggled(bool)));

  connect(mRandomForest, SIGNAL(toggled(bool)), this, SLOT(RandomForestToggled(bool)));

  cvLabel->setEnabled(false);
  mSplitModelDirectoryLabel->setEnabled(false);
  cvValue->setEnabled(false);
  mSplitModelDirectory->setEnabled(false);
  mSplitModelDirectoryButton->setEnabled(false);

  // Disable Optimized ML and FS until they're ready
  mSelectOptimalFS->setVisible(false);
  mSelectOptimalFS->setEnabled(false);
  mSelectOptimalML->setVisible(false);
  mSelectOptimalML->setVisible(false);

  cvValue->setText("5");

  m_jobIsCurrentlyRunning = false;
}
fTrainingSimulator::~fTrainingSimulator()
{
    //m_workerThread->terminate();// Unsafe, but necessary right now, because
    //requestInterruption and quit/exit/wait aren't handled well at the moment
    //m_workerThread->wait();

    //delete m_workerThread;

    //m_workerThread->requestInterruption(); // Advise code to quit ASAP
    //m_workerThread->quit(); // Instruct worker thread to exit
    //m_workerThread->wait(); // Ensure worker thread is actually terminated before continuing
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
void fTrainingSimulator::SVMToggled(bool on)
{
    paramsGroupBox->setVisible(true);
    if (mRandomForestFS->isChecked())
    {
        randomForestGroupBox->setVisible(true);
    }
    else
    {
        randomForestGroupBox->setVisible(false);
    }
}
void fTrainingSimulator::OptimizationToggled(bool on)
{
    if (on)
    {
        paramsGroupBox->setVisible(true);
        randomForestGroupBox->setVisible(false);
        cMinimumSpinbox->setEnabled(true); 
        cMaximumSpinbox->setEnabled(true);
        gMinimumSpinbox->setEnabled(true);
        gMaximumSpinbox->setEnabled(true);
        cMinimumSpinbox->setVisible(true);
        cMaximumSpinbox->setVisible(true);
        gMinimumSpinbox->setVisible(true);
        gMaximumSpinbox->setVisible(true);

        cMinimumLabel->setVisible(true);
        cMaximumLabel->setVisible(true);
        gMinimumLabel->setVisible(true);
        gMaximumLabel->setVisible(true);
    }
    else
    {
        // Just for indicating to the user -- the values still get passed but should remain unused
        cMinimumSpinbox->setEnabled(false);
        cMaximumSpinbox->setEnabled(false);
        gMinimumSpinbox->setEnabled(false);
        gMaximumSpinbox->setEnabled(false);
        cMinimumSpinbox->setVisible(false);
        cMaximumSpinbox->setVisible(false);
        gMinimumSpinbox->setVisible(false);
        gMaximumSpinbox->setVisible(false);

        cMinimumLabel->setVisible(false);
        cMaximumLabel->setVisible(false);
        gMinimumLabel->setVisible(false);
        gMaximumLabel->setVisible(false);
    }
}

void fTrainingSimulator::RandomForestToggled(bool on)
{
    if (on)
    {
        randomForestGroupBox->setVisible(true);
        randomForestEpsilonLabel->setVisible(true);
        randomForestEpsilonSpinbox->setVisible(true);
        randomForestEpsilonSpinbox->setEnabled(true);
        randomForestMaxIterationsLabel->setVisible(true);
        randomForestMaxIterationsSpinbox->setVisible(true);
        randomForestMaxIterationsSpinbox->setEnabled(true);
        paramsGroupBox->setVisible(false);
    }
    else
    {
        if (!mRandomForestFS->isChecked())
        {
            randomForestGroupBox->setVisible(false);
            randomForestEpsilonLabel->setVisible(false);
            randomForestEpsilonSpinbox->setVisible(false);
            randomForestEpsilonSpinbox->setEnabled(false);
            randomForestMaxIterationsLabel->setVisible(false);
            randomForestMaxIterationsSpinbox->setVisible(false);
            randomForestMaxIterationsSpinbox->setEnabled(false);
        }
    }
}
void fTrainingSimulator::CancelButtonPressed()
{
    if (m_jobIsCurrentlyRunning)
    {
        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this, "Confirm cancellation", 
            "Are you sure you want to stop the currently running training job?", QMessageBox::Yes|QMessageBox::No);
        if (reply == QMessageBox::Yes)
        {
            //m_jobIsCurrentlyRunning = false;
            //m_workerThread->terminate();
            //m_workerThread->wait();
            
        }
    }
    else
    {
        this->close();
    }
  

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
  if (mRandomForestFS->isChecked() && mNumberOfFeaturesSpinbox->value() == 0)
  {
      ShowErrorMessage("Please specify a number of features to select. '0' for selecting the best-performing set is currently not supported for Random Forest-based feature selection.", this);
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
    if (!mLinearKernel->isChecked() && !mRBFKernel->isChecked() &&
        !mSigmoidKernel->isChecked() && !mPolynomialKernel->isChecked() &&
        !mChiSquaredKernel->isChecked() && !mIntersectionKernel->isChecked() &&
        !mRandomForest->isChecked() && !mSGDSVM->isChecked() &&
        !mBoostedTrees->isChecked())
    {
      ShowErrorMessage("Please select a classifier.");
      return;
    }
    if (!mSVMFFS->isChecked() && !mESFS->isChecked() &&
        !mBackwardsFS->isChecked() && !mRandomForestFS->isChecked() &&
        !mReliefFFS->isChecked())
    {
      ShowErrorMessage("Please select one of the given feature selection methods.");
      return;
    }
    if ( (mLinearKernel->isChecked() || mRBFKernel->isChecked() ||
        mSigmoidKernel->isChecked() || mPolynomialKernel->isChecked() ||
        mChiSquaredKernel->isChecked() || mIntersectionKernel->isChecked() ) &&
        (mResubstitution->isChecked() == false && mFiveFold->isChecked() == false))
    {
      ShowErrorMessage("Please select at least one of the given two internal cross-validation options: Resubstitution, 5-fold.");
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
    if (mTestPredictionsAgainstProvidedLabels->isChecked() && 
        inputTargetName->text().isEmpty())
     {
        ShowErrorMessage("'Test predictions against given labels' is selected -- Please provide a label file.", this);
        return;
     }


    if (mTestPredictionsAgainstProvidedLabels->isChecked() &&
        !cbica::isFile(inputTargetName->text().toStdString()))
    {
        ShowErrorMessage("The provided label file does not exist -- please double-check.", this);
        return;
    }

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

  // Defaults
  // TODO: Pull these FROM the TrainingModuleParameters class instead of setting them here --  much cleaner
  // Need to update the UI from those too.
  int classifierType = 1;
  int featureselectionType = 1;
  int optimizationType = 0;
  int crossvalidationType = 1;
  int foldType = 10;
  int confType = 1;
  int numberOfFeatures = 0;

  std::string modelpath ="";

  TrainingModuleParameters params; // parameter object passed through

  if (mLinearKernel->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_LINEAR;
  else if (mRBFKernel->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_RBF;
  else if (mPolynomialKernel->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_POLYNOMIAL;
  else if (mSigmoidKernel->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_SIGMOID;
  else if (mIntersectionKernel->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_INTERSECTION;
  else if (mChiSquaredKernel->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SVM_CHISQUARED;
  else if (mRandomForest->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_RANDOMFOREST;
  else if (mSGDSVM->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_SGD_SVM;
  else if (mBoostedTrees->isChecked())
      params.classifierType = CAPTK::ClassifierType::CLASS_TYPE_BOOSTEDTREES;

  if (mSVMFFS->isChecked())
      params.featureSelectionType = CAPTK::FeatureSelectionType::FS_TYPE_FFS;
  else if (mESFS->isChecked())
      params.featureSelectionType = CAPTK::FeatureSelectionType::FS_TYPE_ES;
  else if (mBackwardsFS->isChecked())
      params.featureSelectionType = CAPTK::FeatureSelectionType::FS_TYPE_BACKWARDS;
  else if (mRandomForestFS->isChecked())
      params.featureSelectionType = CAPTK::FeatureSelectionType::FS_TYPE_RANDOMFOREST;
  else if (mReliefFFS->isChecked())
      params.featureSelectionType = CAPTK::FeatureSelectionType::FS_TYPE_RELIEF_F;


  if (mResubstitution->isChecked())
    params.crossValidationType = CAPTK::CrossValidationType::CV_TYPE_RESUBSTITUTION;
  else if (mFiveFold->isChecked())
    params.crossValidationType = CAPTK::CrossValidationType::CV_TYPE_FiveFold;

  if (mOptimization->isChecked())
    params.optimizationType = CAPTK::OptimizationType::OPT_TYPE_ON;
  else
    params.optimizationType = CAPTK::OptimizationType::OPT_TYPE_OFF;

  if (mCrossValidation->isChecked())
  {
    params.configurationType = CAPTK::ClassificationConfigurationType::CONF_TYPE_KFOLD_CV;
    params.folds = cvValue->text().toInt();
  }
  else if (mSplitTrain->isChecked())
  {
    params.configurationType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TRAIN;
  }
  else if (mSplitTest->isChecked())
  {
    params.configurationType = CAPTK::ClassificationConfigurationType::CONF_TYPE_SPLIT_TEST;
  }

  if (mTestPredictionsAgainstProvidedLabels->isChecked())
  {
    params.testPredictionsAgainstProvidedLabels = true;
  }

  params.inputFeaturesFile = mInputFeaturesName.toStdString();
  params.inputLabelsFile = mInputTargetName.toStdString();
  params.outputDirectory = mOutputPathName.toStdString(); 
  //params.modelDirectory = mModelDirectoryName.toStdString();
  params.modelDirectory = mSplitModelDirectory->text().toStdString();
  params.cMin = cMinimumSpinbox->value();
  params.cMax = cMaximumSpinbox->value();
  params.gMin = gMinimumSpinbox->value();
  params.gMax = gMaximumSpinbox->value();
  params.randomForestEpsilon = randomForestEpsilonSpinbox->value();
  params.randomForestMaxIterations = randomForestMaxIterationsSpinbox->value();
  params.maxNumberOfFeatures = mNumberOfFeaturesSpinbox->value();

  
  // TODO fix this to allow threaded in the background with progress updates
  if (!m_jobIsCurrentlyRunning)
  {

      /*
      auto toBeExecuted = std::bind(&fTrainingSimulator::workerFunction, this, params, std::ref(m_lastResult));
      m_workerThread = QThread::create(toBeExecuted);

      connect(m_workerThread, &QThread::finished, m_workerThread, &QObject::deleteLater);
      connect(m_workerThread, SIGNAL(finished()), this, SLOT(onThreadFinished()));

      m_workerThread->start();
      */
      QMessageBox::StandardButton reply;
      reply = QMessageBox::question(this, "Training Module",
          "While the training module is running, you will not be able to interact with the CaPTk interface. \
          Press OK to continue.", QMessageBox::Ok | QMessageBox::Cancel);
      if (reply == QMessageBox::Ok) 
      {
          progressBar->setMinimum(0);
          progressBar->setMaximum(0);
          progressBar->setValue(1);
          progressBar->setTextVisible(true);
          progressBar->setFormat("Running...");
          //ShowMessage("Training Module is now running.", this, "Running...");
          TrainingModule m_trainingSimulator;
          m_jobIsCurrentlyRunning = false;
          m_lastResult.success = false;
          m_lastResult.message = "Operation incomplete!"; // Placeholder for if thread terminates without completion
          m_jobCompleted = false;
          m_lastResult = m_trainingSimulator.Run(params);
          //m_trainingSimulator.RunThread(params);
          m_jobIsCurrentlyRunning = false;
          progressBar->setValue(0);
          progressBar->setMaximum(1);
          if (!m_lastResult.success)
          {
              ShowErrorMessage("Training module failed with the following error: " + m_lastResult.message +
                  "\n Please check parameters and submit a bug report to CBICA Software if applicable.", this, "Error");
          }
          else 
          {
              ShowMessage("Training module finished!", this);
          }
          
      }
      else
      {
          return;
      }

      
  }
  else
  {
      // TODO: Add interruption routine here as needed
      ShowMessage("The training module is already running -- please wait for it to complete.", this, "Please Wait");
  }


}

//void fTrainingSimulator::workerFunction(TrainingModuleParameters& params, TrainingModuleResult& storeResult)
//{
    //storeResult = m_trainingSimulator.Run(params);
    //m_jobCompleted = true;
//}

/*
void fTrainingSimulator::onThreadFinished()
{
    m_jobIsCurrentlyRunning = false;
    progressBar->setMinimum(0);
    progressBar->setMaximum(1);
    progressBar->setValue(0);
    progressBar->reset();
    if (!m_jobCompleted)
    {
        ShowMessage("DEBUG THREAD INTERRUPTED", this);
        // exit silently -- this only is reached if the thread was terminated before exiting TrainingModule::Run()
        logger.Write("fTrainingDialog: User closed the application before the TrainingModule job finished.");
    }
    if (m_lastResult.success)
    {
        QString msg;
        msg = "Trained model has been saved at the specified location.";
        ShowMessage(msg.toStdString(), this);
    }

}
*/

/*
void fTrainingSimulator::onProgressUpdate(int number)
{
    progressBar->setValue(number);
}*/

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