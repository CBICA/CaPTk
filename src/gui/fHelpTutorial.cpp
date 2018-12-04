
#include "fHelpTutorial.h"
//#include "CAPTk.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWebEngineView>
#include <QCheckBox>
#include <QApplication>
#include <QDesktopWidget>
#include "CapTkGUIUtils.h"
#include "cbicaLogging.h"


fHelpTutorial::~fHelpTutorial()
{

}
fHelpTutorial::fHelpTutorial()
{

  // system("open -a ~/3d-nii-visualizer/dist/Theia.app");


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

  // m_webView = new QWebView();
  // NEW CHANGES
  m_webView = new QWebEngineView();

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

  toolbar->addStretch();
  toolbar->addWidget(confirmationCheckBox);
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
