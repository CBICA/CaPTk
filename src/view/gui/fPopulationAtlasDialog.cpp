#include "fPopulationAtlasDialog.h"
//#include "fProgressDialog.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

fPopulationAtlasDialog::fPopulationAtlasDialog()
{
  this->setWindowTitle("Population atlases");
  //TBD read this from common place 
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  //this->setModal(true);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  mode = -1;
  //this->setFixedWidth(400);
  //this->setFixedHeight(300);
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));
  connect(inputdirectoryButton, SIGNAL(clicked()), this, SLOT(OpenInputDirectory()));
  connect(inputlabelButton, SIGNAL(clicked()), this, SLOT(OpenInputLabelFile()));
  connect(inputAtlasButton, SIGNAL(clicked()), this, SLOT(OpenInputAtlasFile()));
  connect(outputdirectoryButton, SIGNAL(clicked()), this, SLOT(OpenOutputDirectory()));
}
fPopulationAtlasDialog::~fPopulationAtlasDialog()
{
}
void fPopulationAtlasDialog::CancelButtonPressed()
{
  this->close();
}
void fPopulationAtlasDialog::ConfirmButtonPressed()
{
  if (inputdirectoryName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the input directory.");
    return;
  }
  if (cbica::directoryExists(inputdirectoryName->text().toStdString()) == false)
  {
    ShowErrorMessage("Input directory does not exist.");
    return;
  }
  //--------------------------------------------------------------------------
  if (outputdirectoryName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the output directory.");
    return;
  }

  if (cbica::directoryExists(outputdirectoryName->text().toStdString()) == false)
  {
    if (cbica::createDirectory(outputdirectoryName->text().toStdString()) == false)
    {
      ShowErrorMessage("Output directory can not be created.");
      return;
    }
  }
  //--------------------------------------------------------------------------
  if (inputlabelName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the label file (.csv).");
    return;
  }
  if (cbica::fileExists(inputlabelName->text().toStdString()) == false)
  {
    ShowErrorMessage("Label file (.csv) does not exist.");
    return;
  }
  //--------------------------------------------------------------------------
  if (inputAtlasName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the atlas file (.nii,.nii.gz).");
    return;
  }
  if (cbica::fileExists(inputAtlasName->text().toStdString()) == false)
  {
    ShowErrorMessage("Atlas file (.nii, .nii.gz) does not exist.");
    return;
  }
  //--------------------------------------------------------------------------
  emit GeneratePopualtionAtlas(inputdirectoryName->text().toStdString(), inputlabelName->text().toStdString(), inputAtlasName->text().toStdString(), outputdirectoryName->text().toStdString());
  //ShowErrorMessage("Atlases saved at the specified location.");
  this->close();
}

void fPopulationAtlasDialog::OpenOutputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull() || directory.isEmpty())
    return;
  else
    outputdirectoryName->setText(directory);
}

void fPopulationAtlasDialog::OpenInputDirectory()
{
  QString directory = getExistingDirectory(this, mInputPathName);
  if (directory.isNull() || directory.isEmpty())
    return;
  else
    inputdirectoryName->setText(directory);
}


void fPopulationAtlasDialog::OpenInputLabelFile()
{
  auto file = getExistingFile(this, mInputPathName, "Labels (*.csv)");

  if (file.isNull() || file.isEmpty())
  {
    return;
  }
  else
  {
    inputlabelName->setText(file);
  }

  //QString extensions = /*CSV_EXT*/".csv";
  //QString inputImage = QFileDialog::getOpenFileName(this, tr("Select File"), mInputPathName, extensions, 0, QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);
  //if (inputImage.isNull())
  //	return;
  //else
  //	inputlabelName->setText(inputImage);
}
void fPopulationAtlasDialog::OpenInputAtlasFile()
{
  auto file = getExistingFile(this, mInputPathName);

  if (file.isNull() || file.isEmpty())
  {
    return;
  }
  else
  {
    inputAtlasName->setText(file);
  }
}
