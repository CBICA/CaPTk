
#include "fHelpTutorial.h"
//#include "CAPTk.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWheelEvent>
#include <QWebEngineView>
#include <QCheckBox>
#include <QApplication>
#include <QPushButton>
#include <QLabel>
#include <QSlider>
#include <QDesktopWidget>
#include "CaPTkGUIUtils.h"
#include "cbicaLogging.h"


fHelpTutorial::~fHelpTutorial()
{

}
fHelpTutorial::fHelpTutorial()
{

  // system("open -a ~/3d-nii-visualizer/dist/Theia.app");

  this->m_zoomValue = 100;
  zoomInButton = new QPushButton();
  zoomInButton->setMaximumSize(30, 30);
  zoomInButton->setToolTip("Zoom In");
  zoomInButton->setText("+");

  zoomOutButton = new QPushButton();
  zoomOutButton->setMaximumSize(30, 30);
  zoomOutButton->setToolTip("Zoom Out");
  zoomOutButton->setText("-");

  zoomSlider = new QSlider(Qt::Orientation::Horizontal);
  zoomSlider->setMaximumWidth(100);
  //! actual range of zoom percent value is from 40% to 400%
  //! slider ranges from 0 - 24
  //! slider value 4 = 100%
  //! incremental step = 15
  zoomSlider->setRange(0, 24);
  zoomSlider->setValue(4);
  zoomSlider->setToolTip("Zoom");

  zoomValueLabel = new QLabel();
  zoomValueLabel->setText("100");

  percentLabel = new QLabel();
  percentLabel->setText("%");


  QVBoxLayout *mainLayout = new QVBoxLayout(this);

  m_dataDir = getCaPTkDataDir();
#if CAPTK_PACKAGE_PROJECT
  m_docDir = cbica::normPath(m_dataDir + "/../share/doc/");
#else
  m_docDir = std::string(PROJECT_SOURCE_DIR) + "../docs/html/";
#endif
  m_helpFileFullPath = m_docDir + "/1_credits.html";

  QHBoxLayout *toolbar = new QHBoxLayout();

  confirmationCheckBox = new QCheckBox();
  confirmationCheckBox->setObjectName("confirmationCheckBox");
  confirmationCheckBox->setText("Never show again");
  confirmationCheckBox->setChecked(false);

  connect(confirmationCheckBox, SIGNAL(toggled(bool)), this, SLOT(on_skipTutorialOnNextRun(bool)));
  connect(zoomInButton, SIGNAL(clicked()), this, SLOT(onZoomInBtnClicked()));
  connect(zoomOutButton, SIGNAL(clicked()), this, SLOT(onZoomOutBtnClicked()));
  connect(zoomSlider, SIGNAL(valueChanged(int)), this, SLOT(onZoomSliderMoved(int)));

  // m_webView = new QWebView();
  // NEW CHANGES
  m_webView = new QWebEngineView();

  //! we want to let the webpage take the max space
    //! this also makes sure that the labels do not get 
    //! preferrence to occupy area thereby resolving issue #548
  m_webView->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  if (!cbica::fileExists(m_helpFileFullPath))
  {
    cbica::Logging(loggerFile, "Unable to start help page, file '" + m_helpFileFullPath + "' not found");
  }
  else
  {
    m_webView->load(QUrl((
#ifndef _WIN32
      "file://" +
#endif
      m_helpFileFullPath).c_str()));
  }

  toolbar->addWidget(confirmationCheckBox);
  toolbar->addStretch();
  toolbar->addWidget(zoomOutButton);
  toolbar->addWidget(zoomSlider);
  toolbar->addWidget(zoomInButton);
  toolbar->addWidget(zoomValueLabel);
  toolbar->addWidget(percentLabel);

  mainLayout->addWidget(m_webView);
  mainLayout->addLayout(toolbar);

  this->setWindowTitle("About CaPTk");
  this->setModal(true);
  QRect rec = QApplication::desktop()->screenGeometry();
  this->setMinimumSize(rec.width() * 0.5, rec.height() * 0.5);
  //this->m_webView->setTextSizeMultiplier(QApplication::desktop()->screen()->logicalDpiX() / 96.0);
  //this->show();

  QCoreApplication::processEvents();

}

void fHelpTutorial::onZoomInBtnClicked()
{
  int zoomStep = 15;
  int newZoomValue = this->m_zoomValue + zoomStep;
  this->SetZoom(newZoomValue);
}

void fHelpTutorial::onZoomOutBtnClicked()
{
  int zoomStep = 15;
  int newZoomValue = this->m_zoomValue - zoomStep;
  this->SetZoom(newZoomValue);
}

void fHelpTutorial::onZoomSliderMoved(int value)
{
  //! rescale: converting the slider value into zoomvalue in range 40 to 400
  int newZoomValue = 40 + (value * 15);
  this->SetZoom(newZoomValue);
}

void fHelpTutorial::wheelEvent(QWheelEvent * event)
{
  if (event->modifiers().testFlag(Qt::ControlModifier))
  {
    QPoint numDegrees = event->angleDelta() / 8;
    if (!numDegrees.isNull())
    {
      int zoomStep = numDegrees.y();
      int newZoomValue = this->m_zoomValue + zoomStep;
      this->SetZoom(newZoomValue);
    }
    event->accept();
  }
}

void fHelpTutorial::SetZoom(int zoomValue)
{
  if (zoomValue >= 40 && zoomValue <= 400)
  {
    if (this->m_zoomValue != zoomValue)
    {
      this->m_webView->setZoomFactor(zoomValue / 100.0);
      this->m_zoomValue = zoomValue;
      this->zoomValueLabel->setText(QString::number(zoomValue));
      this->zoomSlider->blockSignals(true);
      //! rescale: converting the zoomvalue in range 40 to 400 to zoom slider value in range 0-24
      int zoomSliderValue = (zoomValue - 40) / 15;
      this->zoomSlider->setValue(zoomSliderValue);
      this->zoomSlider->blockSignals(false);
    }
  }
}