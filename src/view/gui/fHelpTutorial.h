#pragma once

#include <QDialog>

class QWebEngineView;
class QCheckBox;
class QSlider;
class QLabel;
class QPushButton;

QT_BEGIN_NAMESPACE

class fHelpTutorial : public QDialog
{
  Q_OBJECT

public:
  fHelpTutorial();
  ~fHelpTutorial();

private:

  QCheckBox *confirmationCheckBox;
  QPushButton *zoomInButton;
  QPushButton *zoomOutButton;
  QSlider *zoomSlider;
  QLabel *zoomValueLabel;
  QLabel *percentLabel;
  // QWebView* m_webView;
  // NEW CHANGES
  QWebEngineView * m_webView;
  std::string m_dataDir, m_docDir, m_startPage, m_helpFileFullPath;
  int m_zoomValue;


private slots:
void on_skipTutorialOnNextRun(bool flag)
{
  emit skipTutorialOnNextRun(flag);
}

public slots:
void onZoomInBtnClicked();
void onZoomOutBtnClicked();
void onZoomSliderMoved(int);

protected:
	//! handle ctrl+wheel for zoom in/out
	void wheelEvent(QWheelEvent *event);

	//! zoom to value
	void SetZoom(int zoomValue);

signals:

void skipTutorialOnNextRun(bool flag);


};
