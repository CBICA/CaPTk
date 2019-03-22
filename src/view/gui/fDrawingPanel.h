///////////////////////////////////////////////////////////////////////////////////////
// fDrawingPanel.h
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fDrawingPanel_h_
#define _fDrawingPanel_h_

#include "ui_fDrawingPanel.h"

/**
\class fDrawingPanel

\brief This class controls the elements in the drawing panel of the tab
*/
class fDrawingPanel : public QWidget, private Ui::fDrawingPanel
{
  Q_OBJECT

public:
  //! Constructor
  fDrawingPanel(QWidget * parent = 0);

  //! Destructor
  ~fDrawingPanel() {}
  int getSelectedDrawLabel()
  {
    if (labelSelectorBox != NULL)
    {
      return labelSelectorBox->currentIndex()+1;
    }
    return 1;
  }
  int getSelectedDrawSize()
  {
    if (sizeComboBox != NULL)
    {
      return sizeComboBox->currentIndex();
    }
    return 0;
  }
  void enableShapeButton(QPushButton* button)
  {
    for (int i = 0; i < shapeButtons.size(); i++)
    {
      if (shapeButtons[i] == button)
      {
        shapeButtons[i]->setChecked(true);
      }
      else
      {
        shapeButtons[i]->setChecked(false);
      }
    }
  }

  int getCurrentOpacity() { return m_currentOpacity; };

signals :
  void clearMask(int label=-1);
  void UndoButtonClicked();
  void FillButtonClicked(int);
  void shapesButtonClicked(int mode);
  void CurrentBrushSizeChanged(int);
  void CurrentDrawingLabelChanged(int); // multiLabel related change
  void CurrentMaskOpacityChanged(int); // multiLabel related change
  void helpClicked_Interaction(std::string);
  void sig_ChangeLabelValuesClicked(const std::string, const std::string);


public slots :

  //! Enable voxel based erase functionality
  void shapeEraceButtonFunctionality();

  //! Erase selected label
  void clearSelectedLabelButtonFunctionality();
  
  //! Erase all labels
  void clearAllLabelButtonFunctionality();
  
  //! Change in the size of drawing/erasing brush
  void CurrentSizeChanged(int size);

  //! Change the current selected label
  void CurrentLabelChanged(int size);

  //! Change the current mask opacity
  void CurrentOpacityChanged(int size);

  //! Undo drawing actions
  void UndoButtonFunctionality(); 

  //! Fill button action
  void FillButtonFunctionality();

  //! Draw Shapes
  void shapesNoneButtonFunctionality();
  void shapesLineButtonFunctionality();
  void shapeFreeHandButtonFunctionality();
  void shapesRectangleButtonFunctionality();
  void shapesCircleButtonFunctionality();
  void shapesSphereButtonFunctionality();
  void helpClicked();

  void ChangeLabelValuesClicked();

private:
  int m_currentOpacity = 1;
};


#endif
