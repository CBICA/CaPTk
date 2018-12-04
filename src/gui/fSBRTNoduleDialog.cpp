#include "fSBRTNoduleDialog.h"
#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CapTkGUIUtils.h"

fSBRTNoduleDialog::fSBRTNoduleDialog()
{
  this->setWindowTitle("SBRT Nodule");
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  //connect(maskFileBrowseButton, SIGNAL(clicked()), this, SLOT(OnMaskFileBrowseButtonClicked()));
  //connect(outputDirBrowseButton, SIGNAL(clicked()), this, SLOT(OnOutputDirectoryBrowseButtonClicked()));
  connect(seedImageBrowseButton, SIGNAL(clicked()), this, SLOT(OnSeedImageBrowseButtonClicked()));
  connect(okButton, SIGNAL(clicked()), this, SLOT(OnOKButtonClicked()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(OnCancelButtonClicked()));

}

fSBRTNoduleDialog::~fSBRTNoduleDialog()
{
}

void fSBRTNoduleDialog::OnCancelButtonClicked()
{
  this->close();
}

void fSBRTNoduleDialog::OnSeedImageBrowseButtonClicked()
{
	auto file = getExistingFile(this, mInputPathName);

	if (!file.isEmpty())
	{
		seedImageName->setText(file);
	}
}

//void fSBRTNoduleDialog::OnOutputDirectoryBrowseButtonClicked()
//{
//	QString directory = getExistingDirectory(this, mInputPathName);
//	if (directory.isNull())
//		return;
//	//else
//		//outputDirName->setText(directory);
//}

void fSBRTNoduleDialog::OnOKButtonClicked()
{
	if (!seedImageName->text().isEmpty())
	{
		if (cbica::fileExists(seedImageName->text().toStdString()) == false)
		{
			ShowErrorMessage("seed image does not exist.");
			return;
		}
	}
	this->close();
	emit SBRTNoduleParamReady(seedImageName->text().toStdString(), labelValueSpinBox->value());
}
 