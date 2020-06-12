///////////////////////////////////////////////////////////////////////////////////////
// fAppDownloadDialog.h
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

#ifndef _fAppDownloadDialog_h_
#define _fAppDownloadDialog_h_


//#include "CAPTk.h"
#include "ui_fAppDownloadDialog.h"
#include <QDialog>
#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QUrl>
#include <QProgressDialog>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QMessageBox>

/**
\class fAppDownloadDialog

\brief This class controls the elements in the DICOM converter
*/
class fAppDownloadDialog : public QDialog, private Ui::fAppDownloadDialog
{
  Q_OBJECT

public:
  fAppDownloadDialog();
  ~fAppDownloadDialog();

  // QString m_exe, m_dataDir, m_modelDir; // contains full path and exe name of dcm2nii
  // std::string m_baseModelDir;

  QString downloadPath;
  QString extractPath;
  QString fullPath;
  QString appName;
  std::string downloadLink;

  QUrl url;
  QNetworkAccessManager *manager;
  QNetworkReply *reply;
  QFile *file;
  bool httpRequestAborted;
  qint64 fileSize;
  QString qInputLink;

  void SetPaths(std::string inputPath)
  {
    downloadPath = QString::fromStdString(inputPath);
    // extractPath = QString::fromStdString(inputPath + currentApp);
    extractPath = QString::fromStdString(inputPath);
  }

  void SetDownloadLink(std::string inputLink) {
    qInputLink = QString::fromStdString(inputLink);
    url.setUrl(qInputLink);
  }

  void startRequest(QUrl url);

public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  // slot for readyRead() signal
  void httpReadyRead();

  // slot for finished() signal from reply
  void httpDownloadFinished();

  // slot for downloadProgress()
  void updateDownloadProgress(qint64, qint64);

  void cancelDownload();

signals:
  void doneDownload(QString fullPath, QString extractPath);
};

#endif