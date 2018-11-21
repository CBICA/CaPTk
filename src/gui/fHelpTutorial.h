#pragma once

#include "CAPTk.h"
//#include "ui_fHelpDialog.h"
#include "qlayout.h"
// #include "qwebview.h"
// NEW CHANGES
#include "qwebengineview.h"
#include "qcheckbox.h"

QT_BEGIN_NAMESPACE

class fHelpTutorial : public QDialog
{
  Q_OBJECT

public:
  fHelpTutorial();
  ~fHelpTutorial();

private:

  QCheckBox *confirmationCheckBox;
  // QWebView* m_webView;
  // NEW CHANGES
  QWebEngineView * m_webView;
  std::string m_dataDir, m_docDir, m_startPage, m_helpFileFullPath;


private slots:
void on_skipTutorialOnNextRun(bool flag)
{
  emit skipTutorialOnNextRun(flag);
}

signals:

void skipTutorialOnNextRun(bool flag);


};
