#include "fSBRTAnalysisDialog.h"
#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "cbicaLogging.h"

fSBRTAnalysisDialog::fSBRTAnalysisDialog()
{
  this->setWindowTitle("SBRT Analysis");
  setupUi(this);

  //! hide the result pane. this is now shown via messagebox.
  //! still here in case we need to come back to it
  this->survivalLabel->hide();
  this->survivalValue->hide();
  this->nodalfailureLabel->hide();
  this->nodalFailureValue->hide();
  this->setWindowModality(Qt::NonModal);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  connect(okButton, SIGNAL(clicked()), this, SLOT(OnOKButtonClicked()));
  connect(label, SIGNAL(linkActivated(const QString&)), this, SLOT(OnURLClicked(const QString&)));
  connect(modelDirBtn, SIGNAL(clicked()), this, SLOT(OnSelectModelBtnClicked()));

}

fSBRTAnalysisDialog::~fSBRTAnalysisDialog()
{
}

void fSBRTAnalysisDialog::OnOKButtonClicked()
{
	QString selectedModelDir = modelDirLineEdit->text();
	if (selectedModelDir.isEmpty())
	{
		//model directory is not select ( field empty )
		ShowErrorMessage("Please select model directory.");
	}
	else if (cbica::directoryExists(selectedModelDir.toStdString()))
	{
		//directory is selected and exists
		this->accept();
	}
	else
	{
		//check directory
		ShowErrorMessage("Please check the model directory you selected.");
	}
}

void fSBRTAnalysisDialog::OnURLClicked(const QString&)
{
  QString volumeString;
  volumeString = "You are about to download a model trained on PENN data.\n";
  //volumeString += "Please note that this model was created following certain assumptions \n";
  //volumeString += "(described in the paper below), can be used for research purposes only \n\n\n";
  //volumeString += "H.Akbari, et al. Imaging Surrogates of Infiltration Obtained Via Multi-\n";
  //volumeString += "parametric Imaging Pattern Analysis Predict Subsequent Location of\n";
  //volumeString += "Recurrence of Glioblastoma, Neurosurgery. 78(4) : 572 - 80, 2016.\n";
  //volumeString += "DOI: 10.1227 / NEU.0000000000001202\n";

  QMessageBox *box = new QMessageBox(QMessageBox::Question, "Disclaimer", volumeString, QMessageBox::Ok | QMessageBox::Cancel);
  //box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
  box->setWindowModality(Qt::NonModal);
  QCoreApplication::processEvents();
  if (box->exec() == QMessageBox::Ok)
  {
	  emit DownloadUrl(QUrl(m_trainedModelLink.c_str()));
  }
}

void fSBRTAnalysisDialog::SetSurvivalValue(float val)
{
  this->survivalValue->setText(QString::number(val));
}

void fSBRTAnalysisDialog::SetNodalFailureValue(float val)
{
  this->nodalFailureValue->setText(QString::number(val));
}

void fSBRTAnalysisDialog::OnSelectModelBtnClicked()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull())
    return;
  else
    modelDirLineEdit->setText(directory);

  mInputPathName = directory;
}

void fSBRTAnalysisDialog::OnCancelButtonClicked()
{
  this->reject();
}
