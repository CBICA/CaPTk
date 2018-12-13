#include "fLungFieldDialog.h"
#include "fProgressDialog.h"

fLungFieldDialog::fLungFieldDialog()
{
  this->setWindowTitle("Lung Field");
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  connect(maskFileBrowseButton, SIGNAL(clicked()), this, SLOT(OnMaskFileBrowseButtonClicked()));
  connect(outputDirBrowseButton, SIGNAL(clicked()), this, SLOT(OnOutputDirectoryBrowseButtonClicked()));
  connect(okButton, SIGNAL(clicked()), this, SLOT(OnOKButtonClicked()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(OnCancelButtonClicked()));

}

fLungFieldDialog::~fLungFieldDialog()
{
}

void fLungFieldDialog::OnCancelButtonClicked()
{
  this->close();
}

void fLungFieldDialog::OnMaskFileBrowseButtonClicked()
{
	auto file = getExistingFile(this, mInputPathName);

	if (!file.isEmpty())
	{
		maskFileName->setText(file);
	}
}

void fLungFieldDialog::OnOutputDirectoryBrowseButtonClicked()
{
	QString directory = getExistingDirectory(this, mInputPathName);
	if (directory.isNull())
		return;
	else
		outputDirName->setText(directory);
}

void fLungFieldDialog::OnOKButtonClicked()
{
	if (!maskFileName->text().isEmpty())
	{
		if (cbica::fileExists(maskFileName->text().toStdString()) == false)
		{
			ShowErrorMessage("Mask file does not exist.");
			return;
		}
	}
	this->close();
	emit LungFieldParamReady(maskFileName->text().toStdString(), outputDirName->text().toStdString());
}
 