
#include "fHelpTutorial.h"


fHelpTutorial::~fHelpTutorial()
{

}
fHelpTutorial::fHelpTutorial()
{
  QVBoxLayout *mainLayout = new QVBoxLayout(this);

  if (QDir(QApplication::applicationDirPath() + "/../data/").exists()) // packaged binary
  {
    m_docDir = QApplication::applicationDirPath().toStdString() + "/../share/doc/";
  }
  else if (QDir(QApplication::applicationDirPath() + "/../../data/").exists()) // running from project
  {
    m_docDir = QApplication::applicationDirPath().toStdString() + "/../../docs/tutorial/"; // special provision for tutorial section ONLY
  }
  m_helpFileFullPath = m_docDir + "1_credits.html";

  QHBoxLayout *toolbar = new QHBoxLayout();

  confirmationCheckBox = new QCheckBox();
  confirmationCheckBox->setObjectName("confirmationCheckBox");
  confirmationCheckBox->setText("Never show again");
  confirmationCheckBox->setChecked(false);

  connect(confirmationCheckBox, SIGNAL(toggled(bool)), this, SLOT(on_skipTutorialOnNextRun(bool)));

  m_webView = new QWebView();
  
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
   
  if (!cbica::fileExists(m_helpFileFullPath))
  {
    cbica::Logging(loggerFile, "Unable to start help page, file '" + m_helpFileFullPath + "' not found");
  }
  else
  {
    m_webView->load(QUrl(m_helpFileFullPath.c_str()));
  }

  toolbar->addStretch();
  toolbar->addWidget(confirmationCheckBox);
  mainLayout->addWidget(m_webView);
  mainLayout->addLayout(toolbar);

  this->setWindowTitle("Tutorial");
  this->setModal(true);
  QRect rec = QApplication::desktop()->screenGeometry();
  this->setMinimumSize(rec.width() * 0.5, rec.height() * 0.5);
  //this->show();

  QCoreApplication::processEvents();

}
