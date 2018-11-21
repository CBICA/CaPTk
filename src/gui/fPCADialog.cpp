#include "fPCADialog.h"
#include "fProgressDialog.h"

#include "cbicaITKUtilities.h"

fPCAEstimator::fPCAEstimator()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(outputImageButton, SIGNAL(clicked()), this, SLOT(SelectOutputImage()));
}
fPCAEstimator::~fPCAEstimator()
{
}
void fPCAEstimator::CancelButtonPressed()
{
  this->close();
}
void fPCAEstimator::ConfirmButtonPressed()
{
  if ((inputImageName->text().isEmpty()))
  {
    ShowErrorMessage("Please specify the number of desired principal components.");
    return;
  }
  if (outputImageName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output folder.");
    return;
  }
  emit RunPCAEstimation(inputImageName->text().toInt(), outputImageName->text().toStdString());
  
  this->close();
}

void fPCAEstimator::OpenInputImage()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
	  return;
  else
	  inputImageName->setText(directory);
}

void fPCAEstimator::SelectOutputImage()
{
	QString directory = getExistingDirectory(this, mOutputPathName);
	if (directory.isNull())
		return;
	else
		outputImageName->setText(directory);
}
