///////////////////////////////////////////////////////////////////////////////////////
// fHelpDialog.cxx
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

#include "fHelpDialog.h"
//#include "CAPTk.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWebEngineView>
#include <QPushButton>
#include <QApplication>
#include <QDesktopWidget>
#include "CaPTkGUIUtils.h"
#include "cbicaUtilities.h"
#include "cbicaLogging.h"

fHelpDialog::~fHelpDialog()
{
  m_webView->close();
}
fHelpDialog::fHelpDialog()
{
  QVBoxLayout *mainLayout = new QVBoxLayout(this);

  m_dataDir = getCaPTkDataDir();
#if CAPTK_PACKAGE_PROJECT
  m_docDir = cbica::normPath(m_dataDir + "/../share/doc/");
#else
#if WIN32
  m_docDir = cbica::normPath(captk_currentApplicationPath + "/../docs/html/");
#else
  m_docDir = cbica::normPath(captk_currentApplicationPath + "/../docs/html/");
#endif
#endif
  m_helpFileFullPath = cbica::normPath(m_docDir + "/index.html");

  QHBoxLayout *toolbar = new QHBoxLayout();

  // m_webView = new QWebView();
  // NEW CHANGES
  m_webView = new QWebEngineView();

  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);
  
  homeButton = new QPushButton();
  homeButton->setMaximumSize(30, 30);
  std::string homeIconFullPath = m_dataDir + "icons/home.png";
  QImage homeimage(homeIconFullPath.c_str());
  homeButton->setIcon(QPixmap::fromImage(homeimage));
  homeButton->setIconSize(QSize(homeButton->size().width(), homeButton->size().height()));
  homeButton->setToolTip("Home");

  backButton = new QPushButton();
  backButton->setMaximumSize(30, 30);
  std::string backIconFullPath = m_dataDir + "icons/back.png";
  QImage backimage(backIconFullPath.c_str());
  backButton->setIcon(QPixmap::fromImage(backimage));
  backButton->setIconSize(QSize(homeButton->size().width(), homeButton->size().height()));
  backButton->setToolTip("Back");

  forwardButton = new QPushButton();
  forwardButton->setMaximumSize(30, 30);
  std::string forwardIconFullPath = m_dataDir + "icons/forward.png";
  QImage forwardimage(forwardIconFullPath.c_str());
  forwardButton->setIcon(QPixmap::fromImage(forwardimage));
  forwardButton->setIconSize(QSize(homeButton->size().width(), homeButton->size().height()));
  forwardButton->setToolTip("Forward");

  refreshButton = new QPushButton();
  refreshButton->setMaximumSize(30, 30);
  std::string refreshIconFullPath = m_dataDir + "icons/refresh.png";
  QImage refreshimage(refreshIconFullPath.c_str());
  refreshButton->setIcon(QPixmap::fromImage(refreshimage));
  refreshButton->setIconSize(QSize(homeButton->size().width(), homeButton->size().height()));
  refreshButton->setToolTip("Refresh");

  connect(homeButton,    SIGNAL(clicked()), this, SLOT(on_homeButton_clicked()));
  connect(backButton,    SIGNAL(clicked()), this, SLOT(on_backButton_clicked()));
  connect(forwardButton, SIGNAL(clicked()), this, SLOT(on_forwardButton_clicked()));
  connect(refreshButton, SIGNAL(clicked()), this, SLOT(on_refreshButton_clicked()));

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

  toolbar->addWidget(backButton);
  toolbar->addWidget(forwardButton);
  toolbar->addWidget(refreshButton);
  toolbar->addWidget(homeButton);
  toolbar->addSpacing(800);
  mainLayout->addLayout(toolbar);
  mainLayout->addWidget(m_webView);

  
  this->setWindowTitle("Help");
  //this->setWindowModality(Qt::NonModal); // this needs to be explicitly non-modal to ensure UI stays responsive
  //this->setAttribute(Qt::WA_DeleteOnClose);
  QRect rec = QApplication::desktop()->screenGeometry();
  this->setMinimumSize(rec.width() * 0.5, rec.height() * 0.5);
  //this->m_webView->setTextSizeMultiplier(QApplication::desktop()->screen()->logicalDpiX() / 96.0);

  QCoreApplication::processEvents();

}
void fHelpDialog::on_homeButton_clicked()
{
  m_webView->load(QUrl((
#ifndef _WIN32
      "file://" +
#endif
    m_helpFileFullPath).c_str()));
}
void fHelpDialog::on_backButton_clicked()
{
  m_webView->back();

}
void fHelpDialog::on_forwardButton_clicked()
{
  m_webView->forward();
}

void fHelpDialog::on_refreshButton_clicked()
{
  m_webView->reload();
}

void fHelpDialog::setNewStartPage(const std::string &startPageHTML)
{
  m_dataDir = getCaPTkDataDir();
#if CAPTK_PACKAGE_PROJECT
  m_docDir = cbica::normPath(m_dataDir + "/../share/doc/");
#else
#if _WIN32
  m_docDir = cbica::normPath(captk_currentApplicationPath + "/../docs/html/");
#else
  m_docDir = cbica::normPath(captk_currentApplicationPath + "/docs/html/");
#endif
#endif

  std::string currentStartPage = m_docDir + "/" + startPageHTML;

  if (!cbica::fileExists(m_helpFileFullPath))
  {
    cbica::Logging(loggerFile, "Unable to start help page, file '" + currentStartPage + "' not found");
  }
  else
  {
    m_webView->load(QUrl((
#ifndef _WIN32
      "file://" +
#endif
      currentStartPage).c_str()));
  }

}

