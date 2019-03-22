///////////////////////////////////////////////////////////////////////////////////////
// fDrawingPanel.cpp
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
// Contact details: software@cbica.upenn.edu
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#include "fDrawingPanel.h"
#include "fMainWindow.h"

fDrawingPanel::fDrawingPanel(QWidget * parent) : QWidget(parent)
{
  setupUi(this);

  connect(clearSelectedLabelButton, SIGNAL(clicked()), this, SLOT(clearSelectedLabelButtonFunctionality()));
  connect(clearAllLabelButton, SIGNAL(clicked()), this, SLOT(clearAllLabelButtonFunctionality()));
  connect(sizeComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(CurrentSizeChanged(int)));
  connect(labelSelectorBox, SIGNAL(currentIndexChanged(int)), this, SLOT(CurrentLabelChanged(int)));
  connect(maskOpacitySpinBox, SIGNAL(currentIndexChanged(int)), this, SLOT(CurrentOpacityChanged(int)));
  connect(UndoButton, SIGNAL(clicked()), this, SLOT(UndoButtonFunctionality()));
  connect(shapeFillButton, SIGNAL(clicked()), this, SLOT(FillButtonFunctionality()));
  connect(shapeNoneButton, SIGNAL(clicked()), this, SLOT(shapesNoneButtonFunctionality()));
  connect(shapeEracerButton, SIGNAL(clicked()), this, SLOT(shapeEraceButtonFunctionality()));
  connect(shapeFreeHandButton, SIGNAL(clicked()), this, SLOT(shapeFreeHandButtonFunctionality()));
  connect(shapesLineButton, SIGNAL(clicked()), this, SLOT(shapesLineButtonFunctionality()));
  connect(shapesRectangleButton, SIGNAL(clicked()), this, SLOT(shapesRectangleButtonFunctionality()));
  connect(shapesCircleButton, SIGNAL(clicked()), this, SLOT(shapesCircleButtonFunctionality()));
  connect(shapesSphereButton, SIGNAL(clicked()), this, SLOT(shapesSphereButtonFunctionality()));
  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
  connect(changeButton, SIGNAL(clicked()), this, SLOT(ChangeLabelValuesClicked()));
}

void fDrawingPanel::helpClicked()
{
  emit helpClicked_Interaction("gs_drawing.html");
}

void fDrawingPanel::ChangeLabelValuesClicked()
{
  if (changeOldValues->text().isEmpty())
  {
    ShowErrorMessage("Old values cannot be empty: use format AxBxC with the same number of inputs as new values", this);
    return;
  }
  if (changeNewValues->text().isEmpty())
  {
    ShowErrorMessage("New values cannot be empty: use format AxBxC with the same number of inputs as old values", this);
    return;
  }
  emit sig_ChangeLabelValuesClicked(changeOldValues->text().toStdString(), changeNewValues->text().toStdString());
}

void fDrawingPanel::clearSelectedLabelButtonFunctionality()
{
  emit clearMask(getSelectedDrawLabel());
}
void fDrawingPanel::clearAllLabelButtonFunctionality()
{
  emit clearMask();
}
void fDrawingPanel::CurrentSizeChanged(int size)
{
  emit CurrentBrushSizeChanged(size);
}

void fDrawingPanel::CurrentLabelChanged(int size)
{
  emit CurrentDrawingLabelChanged(size);
}

void fDrawingPanel::CurrentOpacityChanged(int size)
{
  m_currentOpacity = size;
  emit CurrentMaskOpacityChanged(size);
}
void fDrawingPanel::UndoButtonFunctionality()
{
  emit UndoButtonClicked();
}

void fDrawingPanel::shapesNoneButtonFunctionality()
{
  enableShapeButton(shapeNoneButton);
  emit shapesButtonClicked(SHAPE_MODE_NONE);
}
void fDrawingPanel::shapeEraceButtonFunctionality()
{
  enableShapeButton(shapeEracerButton);
  emit shapesButtonClicked(SHAPE_MODE_ERASER);
}
void fDrawingPanel::shapeFreeHandButtonFunctionality()
{
  enableShapeButton(shapeFreeHandButton);
  emit shapesButtonClicked(SHAPE_MODE_FREE_HAND);
}
void fDrawingPanel::shapesLineButtonFunctionality()
{
  enableShapeButton(shapesLineButton);
  emit shapesButtonClicked(SHAPE_MODE_LINE);
}
void fDrawingPanel::shapesRectangleButtonFunctionality()
{
  enableShapeButton(shapesRectangleButton);
  emit shapesButtonClicked(SHAPE_MODE_RECTANGLE);
}
void fDrawingPanel::shapesCircleButtonFunctionality()
{
  enableShapeButton(shapesCircleButton);
  emit shapesButtonClicked(SHAPE_MODE_CIRCLE);
}
void fDrawingPanel::shapesSphereButtonFunctionality()
{
  enableShapeButton(shapesSphereButton);
  emit shapesButtonClicked(SHAPE_MODE_SPHERE);
}

void fDrawingPanel::FillButtonFunctionality()
{
  enableShapeButton(shapeFillButton);
  emit shapesButtonClicked(SHAPE_MODE_FILL);
}
