#include "fPerfusionMeasuresDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

#include "cbicaITKUtilities.h"

fPerfusionEstimator::fPerfusionEstimator()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  m_baselineStart->setValue(0);
  m_baselineEnd->setValue(20);
  m_recoveryStart->setValue(66);
  m_recoveryEnd->setValue(88);
  m_baselineStart->setEnabled(true);
  m_baselineEnd->setEnabled(true);
  m_recoveryStart->setEnabled(true);
  m_recoveryEnd->setEnabled(true);

  m_baselineStartLabel->setText("Baseline Start Threshold %");
  m_baselineEndLabel->setText("Baseline End Threshold %");
  m_recoveryStartLabel->setText("Recovery Start Threshold %");
  m_recoveryEndLabel->setText("Recovery End Threshold %");

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(inputImageButton, SIGNAL(clicked()), this, SLOT(SelectInputImage()));

  m_rcbv->setChecked(true);
  m_psr->setChecked(true);
  m_ph->setChecked(true);
}
fPerfusionEstimator::~fPerfusionEstimator()
{
}
void fPerfusionEstimator::CancelButtonPressed()
{
  this->close();
}
void fPerfusionEstimator::ConfirmButtonPressed()
{
  if ((inputImageName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the DSC-MRI Image.");
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output folder.");
    return;
  }
  if (m_rcbv->isChecked() == false && m_psr->isChecked() == false && m_ph->isChecked() == false)
  {
    ShowErrorMessage("Please select at least one of the given three options: ap-rCBV, PH, PSR.");
    return;
  }
  if (m_baselineStart->value() >= m_baselineEnd->value() || m_recoveryStart->value() >= m_recoveryEnd->value() || m_baselineEnd->value() >= m_recoveryStart->value())
  {
    ShowErrorMessage("Please check your baseline and recovery thresholds for validity.");
    return;
  }
  emit RunPerfusionMeasuresCalculation(m_rcbv->isChecked(), m_psr->isChecked(), m_ph->isChecked(), m_baselineStart->value(), m_baselineEnd->value(), m_recoveryStart->value(), m_recoveryEnd->value(), mInputPathName.toStdString(), mOutputPathName.toStdString());

  this->close();
}

void fPerfusionEstimator::SelectOutputImage()
{
  QString directory = getExistingDirectory(this, mOutputPathName);
  if (directory.isNull() || directory.isEmpty())
    return;
  else
    outputImageName->setText(directory);

  mOutputPathName = directory;
}

void fPerfusionEstimator::SelectInputImage()
{
	auto inputImage = getExistingFile(this, mInputPathName);
	if (inputImage.isNull() || inputImage.isEmpty())
		return;
	else
		inputImageName->setText(inputImage);

	mInputPathName = inputImage;
}
