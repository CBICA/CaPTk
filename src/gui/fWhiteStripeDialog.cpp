#include "fWhiteStripeDialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"
#include "CaPTkGUIUtils.h"

fWhiteStripeObj::fWhiteStripeObj()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;

  //options_T1selected->setChecked(true);
  //options_T2selected->setChecked(false);
  options_skullStrippedImage->setChecked(true);
  options_axialSlicing->setChecked(false);
  outputImageName->setText(mInputPathName + "/whiteStripe_output.nii.gz");

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
  connect(options_T1Selector, SIGNAL(currentIndexChanged(int)), this, SLOT(T1ComboBoxInvoked(int)));
  //connect(options_T1selected, SIGNAL(toggled(bool)), this, SLOT(T1ImageSelected()));
  //connect(options_T2selected, SIGNAL(toggled(bool)), this, SLOT(T2ImageSelected()));
  connect(options_skullStrippedImage, SIGNAL(toggled(bool)), this, SLOT(SetSkullStrippedImage()));
  connect(options_axialSlicing, SIGNAL(toggled(bool)), this, SLOT(SetAxialSlicingNeeded()));
}
fWhiteStripeObj::~fWhiteStripeObj()
{
}
void fWhiteStripeObj::CancelButtonPressed()
{
  this->close();
}
void fWhiteStripeObj::ConfirmButtonPressed()
{
  if (outputImageName->text().isEmpty())
  {
    outputImageName->setText(mInputPathName + "whiteStripe_output.nii.gz");
  }

  //bool T1selected = true;

  //if (options_T2selected->isEnabled())
  //{
  //  T1selected = false;
  //}

  emit RunWhiteStripe(options_twsWidth->text().toDouble(), options_sliceStartZ->text().toInt(), options_sliceStopZ->text().toInt(),
    options_tissuesMax->text().toInt(), options_smoothMax->text().toDouble(), options_smoothDelta->text().toDouble(),
    options_histSize->text().toInt(), t1Image, outputImageName->text().toStdString());
  
  this->close();
}

void fWhiteStripeObj::T1ComboBoxInvoked(int index)
{
  if (options_T1Selector->currentIndex() == 0)
  {
    t1Image = true;
  }
  else
  {
    t1Image = false;
  }
}

//void fWhiteStripeObj::T1ImageSelected()
//{
//  options_T2selected->setChecked(false);
//  t1Image = true;
//  //checkSkullStripOption();
//}
//
//void fWhiteStripeObj::T2ImageSelected()
//{
//  options_T1selected->setChecked(false);
//  t1Image = false;
//  //checkSkullStripOption();
//}

void fWhiteStripeObj::SetSkullStrippedImage()
{
  options_axialSlicing->setChecked(false);
  skullStrippedImage = true;
  options_sliceStartZ->setDisabled(true);
  options_sliceStopZ->setDisabled(true);
  options_sliceStartZ->setText("-1");
  options_sliceStopZ->setText("-1");
  //checkT1SelectedOption();
}

void fWhiteStripeObj::SetAxialSlicingNeeded()
{
  options_skullStrippedImage->setChecked(false);
  skullStrippedImage = false;
  options_sliceStartZ->setDisabled(false);
  options_sliceStopZ->setDisabled(false);
  options_sliceStartZ->setText("80");
  options_sliceStopZ->setText("120");
  //checkT1SelectedOption();
}

void fWhiteStripeObj::SelectOutputImage()
{
  QString outputImage = getSaveFile(this, mInputPathName, mInputPathName + "whiteStripe_output.nii.gz");
  if (outputImage.isNull())
    return;
  else
    outputImageName->setText(outputImage);

  mInputPathName = cbica::getFilenameBase(outputImage.toStdString(), false).c_str(); // overwrite previous default path with new output
}
