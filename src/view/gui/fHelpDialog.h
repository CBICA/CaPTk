///////////////////////////////////////////////////////////////////////////////////////
// fHelpDialog.h
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

#ifndef _fHelpDialog_h_
#define _fHelpDialog_h_

#include <QDialog>

class QWebEngineView;
class QSlider;
class QLabel;

QT_BEGIN_NAMESPACE

/**
\class fHelpDialog

\brief This class controls the elements in the help dialog
*/
class fHelpDialog : public QDialog
{
  Q_OBJECT
private:
  // QWebView* m_webView;
  // NEW CHANGES 
  QWebEngineView* m_webView;
  QPushButton *homeButton;
  QPushButton *backButton;
  QPushButton *forwardButton;
  QPushButton *refreshButton;
  QPushButton *zoomInButton;
  QPushButton *zoomOutButton;
  QSlider *zoomSlider;
  QLabel *zoomValueLabel;
  QLabel *percentLabel;
  std::string m_helpFileFullPath;
  std::string m_dataDir, m_docDir;
  std::string m_startPage;
  int zoomPercentValue;
   
private slots:
  void on_homeButton_clicked();
  void on_backButton_clicked();
  void on_forwardButton_clicked();
  void on_refreshButton_clicked();
  void on_zoomInButton_clicked();
  void on_zoomOutButton_clicked();
  void on_zoomSliderMoved(int);

protected:
  //! handle ctrl+wheel for zoom in/out
  void wheelEvent(QWheelEvent *event);

  //! zoom to value
  void SetZoom(int zoomValue);

public:

  fHelpDialog();
  ~fHelpDialog();

  void setNewStartPage(const std::string &startPageHTML);
};
#endif
