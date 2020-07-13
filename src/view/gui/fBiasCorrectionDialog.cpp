#include "fBiasCorrectionDialog.h"
#include "CaPTkGUIUtils.h"
#include "BiasCorrection.hpp"

fBiasCorrectionDialog::fBiasCorrectionDialog()
{
  setupUi(this);
  this->setModal(true);
  confirmButton->setEnabled(false); // wait for user selection

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(mode_radioButtons, SIGNAL(buttonClicked(int)), this, SLOT(SelectedMode(int)));

  this->LoadDefaultParameters();
}

fBiasCorrectionDialog::~fBiasCorrectionDialog()
{
}

void fBiasCorrectionDialog::ConfirmButtonPressed()
{
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output file.", this);
    return;
  }

  // Read text as numbers, check if OK
  bool splineOrder_ok, otsuBins_ok, maxIterations_ok, fittingLevels_ok, filterNoise_ok, fwhm_ok;
  int splineOrder = options_splineOrderName->text().toInt(&splineOrder_ok);
  int otsuBins = options_otsuBinsName->text().toInt(&otsuBins_ok);
  int maxIterations = options_maxIterationsName->text().toInt(&maxIterations_ok);
  int fittingLevels = options_fittingLevelsName->text().toInt(&fittingLevels_ok);
  float filterNoise = options_filterNoiseName->text().toFloat(&filterNoise_ok);
  float fwhm = options_fwhmName->text().toFloat(&fwhm_ok);

  if (!splineOrder_ok || !otsuBins_ok || !maxIterations_ok || !fittingLevels_ok || !filterNoise_ok || !fwhm_ok)
  {
    std::string msg = "Check your parameters: parameters should be numeric. "
      "Spline order, Otsu bins, max iterations, and fitting levels should be integers. "
      "Filter noise and fwhm should be floating-point values.";
    ShowErrorMessage(msg, this);
    return;
  }

  emit CallBiasCorrection(this->correctionMode,
    outputImageName->text(),
    splineOrder,
    otsuBins,
    maxIterations,
    fittingLevels,
    filterNoise,
    fwhm
  );

  this->close();
}

void fBiasCorrectionDialog::CancelButtonPressed()
{
  this->close();
}
void fBiasCorrectionDialog::SelectOutputImage()
{
  QString saveFileName = getSaveFile(this, mInputPathName, mInputPathName + "biasCorrect.nii.gz");
  outputImageName->setText(saveFileName);
}
// if either is selected, enable confirm for the rest of the session
// (impossible to "un-choose" and the setting persists)
void fBiasCorrectionDialog::SelectedMode(int mode)
{
  if (mode == 0)
  {
    this->correctionMode = "N3";
  }
  else if (mode == 1)
  {
    this->correctionMode = "N4";
  }

  this->confirmButton->setEnabled(true);
}

void fBiasCorrectionDialog::LoadDefaultParameters()
{
  this->options_splineOrderName->setText(QString::number(BiasCorrection::default_splineOrder));
  this->options_otsuBinsName->setText(QString::number(BiasCorrection::default_otsuBins));
  this->options_maxIterationsName->setText(QString::number(BiasCorrection::default_maxIterations));
  this->options_fittingLevelsName->setText(QString::number(BiasCorrection::default_fittingLevels));
  this->options_filterNoiseName->setText(QString::number(BiasCorrection::default_filterNoise));
  this->options_fwhmName->setText(QString::number(BiasCorrection::default_fwhm));
}