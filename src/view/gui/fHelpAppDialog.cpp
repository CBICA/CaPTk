///////////////////////////////////////////////////////////////////////////////////////
// fHelpAppDialog.cxx
//
// Copyright (c) 2016. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/cbica/captk/license.html
///////////////////////////////////////////////////////////////////////////////////////



#include "fHelpAppDialog.h"


fHelpAppDialog::fHelpAppDialog()
{
  this->setFixedSize(520, 600);
  setupUi(this);
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  //this->setAttribute(Qt::WA_DeleteOnClose);
  SetLoadingTab();
  SetTumorTab();
  SetShortcutsTab();
  SetDrawingTab();
  this->setWindowModality(Qt::NonModal);
  this->show();
  QCoreApplication::processEvents();
}
fHelpAppDialog::~fHelpAppDialog()
{

}
void fHelpAppDialog::SetLoadingTab()
{

  label1->setText(QString::fromStdString("\nTest_1:"));
  label2->setText(QString::fromStdString("\nUnloading:"));
  label3->setText(QString::fromStdString("\nSaving:"));
  label4->setText(QString::fromStdString("\nOverlay:"));
  label5->setText(QString::fromStdString("\nColor mappings:"));


  label1text->setText(QString::fromStdString("\nMRI (NIFTI): Load all co-registered images using \"File -> Load -> \nMRI (NIFTI)\".\nMRI (NIFTI): Load DICOM images using \"File -> Load -> MRI (DICOM)\". \nROI (NIFTI): Load the existing ROI using \"File -> Load -> ROI (NIFTI)\"."));
  label2text->setText(QString::fromStdString("\nUnload individual images by clicking the \"X\" button on the left\nof each image filename."));
  label3text->setText(QString::fromStdString("\nMRI (NIFTI): Save individual images by using \"File -> Save -> \nMRI (NIFTI)\".\nROI (NIFTI): Save ROI by using \"File -> Save -> ROI (NIFTI)\"."));
  label4text->setText(QString::fromStdString("\nThe \"Overlay\" tickbox  allows for observation of intensity\nchanges from one image to another using a slider. Once this is\nenabled, the user must select one image from the top list and\nanother image from the bottom list of images, which will\ncorrespond to the images on the left and the right side of the\nslider, respectively."));
  label5text->setText(QString::fromStdString("\nThe \"Preset\" and the rest of the controls in the bottom left side\nof the panel, change the color mapping of the visualization\nwindows, in order to depict different properties of the displayed\nimages."));

}

void fHelpAppDialog::SetTumorTab()
{

  label2t->setText(QString::fromStdString("\nTest_2"));
  label3t->setText(QString::fromStdString(""));
  label4t->setText(QString::fromStdString(""));
  label5t->setText(QString::fromStdString("\nLists:"));
  label6t->setText(QString::fromStdString(""));
  label7t->setText(QString::fromStdString(""));




  label2ttext->setText(QString::fromStdString("\nChoose the type of points you want to initialise."));
  label3ttext->setText(QString::fromStdString("\"Tumor points\" are used to approximate the bulk volume of each           .\napparent tumor by a sphere, and they refer to a seed point for a\ntumor center and another point to define the radius of the sphere."));
  label4ttext->setText(QString::fromStdString("\n\"Tissue points\" are used to model the intensity distribution of each\n brain tissue type, e.g. white matter (WM), grey matter (GM), \ncerebellum (CB)."));
  label5ttext->setText(QString::fromStdString("\nThe RAS coordinates of the related points are shown in these lists."));
  label6ttext->setText(QString::fromStdString("Erase individual points by selecting (double-clicking) them in the\nlist and using the \"Remove\" button, or erase all of them by clicking\n the \"Clear all\" button."));
  label7ttext->setText(QString::fromStdString("\n\"Load\"/\"Save\" buttons are used to manipulate the related points \nfrom/to a file. Pre-specified filenames are suggested depending on \nthe selected target application."));


}


void fHelpAppDialog::SetDrawingTab()
{


  drawing->setText(QString::fromStdString("\nTest_3"));
  erasing->setText(QString::fromStdString("\nCLEAR CONTROLS:"));

  label1d->setText(QString::fromStdString("\nNear ROI:"));
  label2d->setText(QString::fromStdString("\nFar ROI:"));
  label3d->setText(QString::fromStdString("\nIndividual voxels:"));
  label4d->setText(QString::fromStdString("\nMarker size:"));

  label5d->setText(QString::fromStdString("\nNear ROI:"));
  label6d->setText(QString::fromStdString("\nFar ROI:"));

  label1dtext->setText(QString::fromStdString("\nTurns the system to drawing mode if it is in view \nmode, and enables the drawing of the Near ROI."));
  label2dtext->setText(QString::fromStdString("\nTurns the system to drawing mode if it is in view \nmode, and enables the drawing of the Near ROI."));
  label3dtext->setText(QString::fromStdString("\nTurns the system to erasing mode, allowing the user \nto erase individual voxels from the Near and/or the \nFar ROI."));
  label4dtext->setText(QString::fromStdString("\nAllows the user to select the size of the marker used \nboth for Drawing or Erasing individual voxels from \nthe Near/Far ROIs."));

  label5dtext->setText(QString::fromStdString("\nClears the Near ROI from all slices."));
  label6dtext->setText(QString::fromStdString("\nClears the Far ROI from all slices."));

}


void fHelpAppDialog::SetShortcutsTab()
{

  descriptionm->setText(QString::fromStdString("Test_4"));
  keyboardheading->setText(QString::fromStdString("\nKEYBOARD:"));
  mouseheading->setText(QString::fromStdString("\nMOUSE:"));

  label1s->setText(QString::fromStdString("\nViewing use:"));
  label2s->setText(QString::fromStdString("\nGeometrical adjustments:"));
  label3s->setText(QString::fromStdString("\nIntensity adjustments:"));

  label4s->setText(QString::fromStdString("\nViewing use:"));
  label5s->setText(QString::fromStdString("\nDrawing mode:"));
  label6s->setText(QString::fromStdString("\nTumor point initialization:"));
  label7s->setText(QString::fromStdString("\nTissue point initialization:"));

  label1stext->setText(QString::fromStdString("\n<U/D>: Change slice number\n<LB> + <movement>: Move crosshair"));
  label2stext->setText(QString::fromStdString("\nCtrl + <U/D>: Zoom in/out (Adjust scaling factor)\nCtrl + <LB> + <movement>: Pan images (Adjust image's center)"));
  label3stext->setText(QString::fromStdString("\n<RB> + <vertical movement>: Brightness adjustment\n<RB> + <horizontal movement>: Contrast adjustment"));

  label4stext->setText(QString::fromStdString("\n1-4: Change view among loaded images\nh: Hide/Display crosshair & seeded points\nr: Reset \"Geometrical Adjustments\"\na: Reset \"Intensity Adjustments\""));
  label5stext->setText(QString::fromStdString("\nf: Enable/Disable foreground drawing\nb: Enable/Disable background drawing\ne: Enable/Disable erase moden\n1-9: Change brush size\nEsc: Quit drawing mode"));
  label6stext->setText(QString::fromStdString("\nShift + spacebar: Set the tumor centre\nCtrl + spacebar: Update the tumor radius\nspacebar: Update the tumor center location"));
  label7stext->setText(QString::fromStdString("\nShift + spacebar: Set a new point for a specified tissue\nCtrl + spacebar: Set a new point for a specified tissue\nspacebar: Update the location of a specified tissue point"));

}