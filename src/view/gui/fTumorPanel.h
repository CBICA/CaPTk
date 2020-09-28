///////////////////////////////////////////////////////////////////////////////////////
// fTumorPanel.h
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
// License Agreement: https://www.med.upenn.edu/cbica/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fTumorPanel_h_
#define _fTumorPanel_h_

//#include "CAPTk.h"
#include "Landmarks.h"
#include "ui_fTumorPanel.h"


/**
\class fTumorPanel

\brief This class controls the elements in the seed initialization dialog of the tab UI
*/
class fTumorPanel : public QWidget, private Ui::fTumorPanel
{
  Q_OBJECT

public:
  //! Default Constructor
  fTumorPanel(QWidget * parent = 0);

  //! Default Destructor
  ~fTumorPanel() {}

  /**
  \brief Start initialization of "seed" points (tumor points)
  */
  void SetCurrentSPoints(Landmarks *lm);

  /**
  \brief Start initialization of tissue points
  */
  void SetCurrentTPoints(Landmarks *lm);

  /**
  \brief Set landmark type (either "seed" or tissue points)
  */
  void SetCurrentSelectedTissueType();

  /**
  \brief Set the current image(s) path
  */
  void SetCurrentPath(const std::string &path)
  {
    mCurrentPath = path;
  }

  /**
  \brief Clear all landmark points
  */
  void Clear();

  /**
  \brief Move slicer cursor to the coordinates pointed in the table for Tissue Table
  */
  void HandleKeyPressingEventTTable();

  /**
  \brief Scroll down in the table when down arrow is pressed for Tissue Table
  */
  void HandleDownKeyEventTTable();

  /**
  \brief Scroll up in the table when up arrow is pressed for Tissue Table
  */
  void HandleUpKeyEventTTable();

  /**
  \brief Delete entry when delete key is pressed for Tissue Table
  */
  void HandleDeleteKeyEventTTable();

  /**
  \brief Move slicer cursor to the coordinates pointed in the table for Seed Table
  */
  void HandleKeyPressingEventSTable();

  /**
  \brief Scroll down in the table when down arrow is pressed for Seed Table
  */
  void HandleDownKeyEventSTable();

  /**
  \brief Scroll up in the table when up arrow is pressed for Seed Table
  */
  void HandleUpKeyEventSTable();

  /**
  \brief Delete entry when delete key is pressed for Seed Table
  */
  void HandleDeleteKeyEventSTable();


  public slots:
  //! gets this signal when panlel tab (Images ,  Draw , Seed etc get changed)
  void tabSelected()
  {
    m_typeRadBtnAllTissues->setChecked(true);
  }
  //! Highlight the selected points in the table
  void HighlightCurrentSelctedPoints(double x, double y, double z, double X, double Y, double Z, double value);

  //! Load seed point(s)
  void sLoad();
  void sLoad(QString file);

  //! Save seed point(s)
  void sSave();

  //! Load tissue point(s)
  void tLoad();
  void tLoad(QString file);

  //! Save tissue point(s)
  void tSave();
  void tSave(QString file);

  //! Remove selected seed point(s)
  void sRemoveSelectedPoints();

  //! Remove selected tissue point(s)
  void tRemoveSelectedPoints();

  //! Remove current tissue point pointed in the row
  void tRemoveCurrentIndexPoint(unsigned int rowIndex);

  //! Remove current seed point pointed in the row
  void sRemoveCurrentIndexPoint(unsigned int rowIndex);

  //! This addition is a wrapper
  void sAddPoint();

  //! Add seed point based on the landmark index and update the table
  void sAddPoint(int landmarksIndex, bool update);
  

  //! Add tissue point based on the landmark index
  void tAddPoint(int rowIndex);

  //! Fix cursor to double clicked entry in seed point table
  void sTableDoubleClicked(int row, int col);

  //! Fix cursor to double clicked entry in tissue point table
  void tTableDoubleClicked(int row, int col);

  //! Highlight the text if seed point table is selected
  void sTableFocused(bool bFocused);

  //! Highlight the text if tissue point table is selected
  void tTableFocused(bool bFocused);

  //! Remove all tissue points
  void tRemoveAllPoints();

  //! Remove all seed points
  void sRemoveAllPoints();

  //! Set seed point type
  void SetSeedType();

  //! Set tissue point type while passing tissue type
  void SetTissueType(int inputTissueType);

  //! Set tissue point type
  void SetTissueType();

  //! Set GLISTR tissue point type
  void SetGLISTRTissueType(); // genericLabels - replace application-specific saving by a single enum-controlled saving mechanism which gets passed to SetTissueType(int)

  //! Set PORTR pre tissue point type
  void SetPORTRPRETissueType(); // genericLabels - replace application-specific saving by a single enum-controlled saving mechanism which gets passed to SetTissueType(int

  //! Set PORTR post tissue point type
  void SetPORTRPOSTTissueType(); // genericLabels - replace application-specific saving by a single enum-controlled saving mechanism which gets passed to SetTissueType(int

  //! Update current points
  void UpdateCurrentPoints();

  void helpClicked();

signals:
  void UpdateRenderWindows();
  void SetActiveLandmarksTypeSignal(int, int, int);
  void MoveSlicerCursor(double, double, double);
  void UpdatePointsRequestByTumorPanel(int, VectorDouble);
  void SetTissueCounter(int type);
  void SetTissueTableIndex(int row);
  void helpClicked_Interaction(std::string);

public:
  int mType;
  int mTissueType;
  Landmarks* mCurrentSPoints;
  Landmarks* mCurrentTPoints;
  std::string mCurrentPath;

  bool mTumorPointsSelected = false;

private:
  std::string initializationFileName; //! file name to save the landmark points

  //! flags for different save types
  enum SaveType
  {
    allTissues, GLISTR, PORTRPre, PORTRPost, Generic
  };

  int m_saveType = SaveType::allTissues;

  //! Populate labels for label radio buttons
  inline void RadioButtonLabels();
};


#endif