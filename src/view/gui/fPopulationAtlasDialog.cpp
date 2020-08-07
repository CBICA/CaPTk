#include "fPopulationAtlasDialog.h"
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
  connect(inputfileButton, SIGNAL(clicked()), this, SLOT(OpenInputFile()));
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
  if (inputfileName->text().isEmpty())
  {
    ShowErrorMessage("Please specify the label file (.csv).");
    return;
  }
  if (cbica::fileExists(inputfileName->text().toStdString()) == false)
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
  emit GeneratePopualtionAtlas(inputfileName->text().toStdString(), inputAtlasName->text().toStdString(), outputdirectoryName->text().toStdString());
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

void fPopulationAtlasDialog::OpenInputFile()
{
  auto file = getExistingFile(this, mInputPathName, "Labels (*.csv)");

  if (file.isNull() || file.isEmpty())
  {
    return;
  }
  else
  {
    inputfileName->setText(file);
  }
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
