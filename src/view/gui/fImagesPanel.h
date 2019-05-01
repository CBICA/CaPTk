///////////////////////////////////////////////////////////////////////////////////////
// fImagesPanel.h
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

#ifndef _fImagesPanel_h_
#define _fImagesPanel_h_

//#include "CAPTk.h"
#include "Landmarks.h"
#include "ui_fImagesPanel.h"
#include "QTablePushButton.h"

enum TABLE_IMAGES_COLUMNS
{
  IMAGES_COLUMN_CLOSE, IMAGES_COLUMN_TYPE, IMAGES_COLUMN_NAME, IMAGES_COLUMN_MODALITY, IMAGES_COLUMN_OVERLAY
};

/**
\class fImagesPanel

\brief This class controls the elements in the images panel of the tab
*/
class fImagesPanel : public QWidget, private Ui::fImagesPanel
{
  Q_OBJECT

public:
  //! Default Constructor
  fImagesPanel(QWidget * parent = 0);

  //! animate comparison mode button click
  void CompareButtonClick();

  //! Default Destructor
  ~fImagesPanel() {}
  void NewImageLoaded(QString idstr, const std::string &filename, int rowIndex, const std::string &imageSubTypeStr, const int imgSubtype, const QObject* caller);
  QTableWidget * GetImagesTable()
  {
    return m_imagesTable;
  }
  QTableWidget * GetNonViewingImagesTable()
  {
    return m_nonVisImagesTable;
  }
  int getModality(const int rowIndex)
  {
    int modality = 0;
    if (m_imagesTable != NULL && rowIndex< m_imagesTable->rowCount())
    {
      QComboBox* cmbModality = (QComboBox*)m_imagesTable->cellWidget(rowIndex, IMAGES_COLUMN_MODALITY);
      if (cmbModality != NULL)
      {
        return cmbModality->currentIndex();
      }
    }
    return modality;
  }
private:


public slots:
  void ImageTableSelectionChanged(QTableWidgetItem*);
  void ImageTableSelectionChanged();
  void overlayUseStateChanged(int);
  void overlaySliderChanged(int);
  void overlayChanged();
  void imageModalityChanged(int);
  void helpClicked();
  void theiaClicked();

  QTableWidgetItem* getSelectedImage()
  {
    QTableWidgetItem* selected = NULL;
    if (m_imagesTable != NULL)
    {
      QList<QTableWidgetItem*> items = m_imagesTable->selectedItems();
      if (!items.empty())
      {
        selected = items[0];
      }
    }
    return selected;
  }
  QTableWidgetItem* getSelectedOverlay()
  {
    if (m_imagesTable == NULL) return NULL;

    QTableWidgetItem* selected = NULL;
    for (int i = 0; i < m_imagesTable->rowCount(); i++)
    {
      QRadioButton *rb = qobject_cast<QRadioButton*>(m_imagesTable->cellWidget(i, 4));
      if (rb == NULL)
      {
        return NULL;
        //TBD fix this issue 
      }
      if (rb->isChecked())
      {
        selected = m_imagesTable->item(i, 2);//TBD check This
      }
    }
    return selected;
  }

signals:
  void sigImageTableSelectionChanged();
  void sigOverlayCheckBoxChanged(int);
  void sigOverlaySliderChanged(int);
  void sigOverlayChanged();
  void sigImageModalityChanged(int);
  void helpClicked_Interaction(std::string);
  void sigTheiaClicked();
  void CompareModeToggled(bool);

};


#endif
