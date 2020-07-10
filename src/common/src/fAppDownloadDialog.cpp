#include <QDialog>

//#include <CaPTkDefines.h>
//#include "CaPTkUtils.h"

#include "fAppDownloadDialog.h"
#include "ApplicationPreferences.h"

fAppDownloadDialog::fAppDownloadDialog()
{
  setupUi(this);
  this->setWindowModality(Qt::NonModal);
  this->setModal(true); // this is a pre-processing routine and therefore should be modal
  this->setWindowFlags(this->windowFlags() & ~Qt::WindowContextHelpButtonHint);

  connect(cancelButton, SIGNAL(clicked()), this, SLOT(CancelButtonPressed()));
  connect(confirmButton, SIGNAL(clicked()), this, SLOT(ConfirmButtonPressed()));

}
fAppDownloadDialog::~fAppDownloadDialog()
{
}

void fAppDownloadDialog::CancelButtonPressed()
{
    this->close();
}

void fAppDownloadDialog::ConfirmButtonPressed()
{
    initDownload();
}

void fAppDownloadDialog::initDownload() {
    ApplicationPreferences::GetInstance()->DeSerializePreferences();
    bool downloadFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadFinishedStatus()).toBool();

    if (downloadFinished) { // if download is done but file is not found and extraction progress is not picked up then reset
        ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SerializePreferences();
    }

    ApplicationPreferences::GetInstance()->DeSerializePreferences();
    bool downloadStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadStartedStatus()).toBool();
    downloadFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadFinishedStatus()).toBool();
    bool extractionStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionStartedStatus()).toBool();
    bool extractionFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionFinishedStatus()).toBool();

    // ApplicationPreferences::GetInstance()->DisplayPreferences();

    if (!downloadStarted)
    {
        setupDownload(this);

        // connect(downloadProgressDialog, SIGNAL(canceled()), this, SLOT(cancelDownload()));

        manager = new QNetworkAccessManager(this);

        // url = 
        QFileInfo fileInfo(url.path());
        QString fileName = fileInfo.fileName();
        appName = fileName;

        if (fileName.isEmpty())
            fileName = "index.html";

        fullPath = downloadPath + fileName;
        // ShowErrorMessage(fullPath.toStdString());

        if (QFile::exists(fullPath)) {
            if (QMessageBox::question(this, tr("HTTP"),
                    tr("There already exists a file %1. Overwrite?").arg(fileName),
                    QMessageBox::Yes|QMessageBox::No, QMessageBox::No)
                    == QMessageBox::No)
                    return;
            QFile::remove(fullPath);
        }

        file = new QFile(fullPath);
        if (!file->open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, tr("HTTP"),
                        tr("Unable to save the file %1: %2")
                        .arg(fileName).arg(file->errorString()));
            delete file;
            file = nullptr;
            return;
        }
        
        // used for progressDialog
        // This will be set true when canceled from progress dialog
        httpRequestAborted = false;

        ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("true").toString());
        ApplicationPreferences::GetInstance()->SerializePreferences();
        // ApplicationPreferences::GetInstance()->DisplayPreferences();
        // downloadProgressDialog->setWindowTitle(tr("HTTP"));
        // downloadProgressDialog->setLabelText(tr("Downloading %1").arg(fileName));
        
        // emit startDownload();

        startRequest(url);

        this->close();
    }
}

// This will be called when download button is clicked
void fAppDownloadDialog::startRequest(QUrl url)
{
    // get() method posts a request
    // to obtain the contents of the target request
    // and returns a new QNetworkReply object
    // opened for reading which emits
    // the readyRead() signal whenever new data arrives.
    reply = manager->get(QNetworkRequest(url));

    // Whenever more data is received from the network,
    // this readyRead() signal is emitted
    connect(reply, SIGNAL(readyRead()),
            this, SLOT(httpReadyRead()));

    // Also, downloadProgress() signal is emitted when data is received
    connect(reply, SIGNAL(downloadProgress(qint64,qint64)),
            this, SLOT(updateDownloadProgress(qint64,qint64)));

    // This signal is emitted when the reply has finished processing.
    // After this signal is emitted,
    // there will be no more updates to the reply's data or metadata.
    connect(reply, SIGNAL(finished()),
            this, SLOT(httpDownloadFinished()));
}


void fAppDownloadDialog::httpReadyRead()
{
    // this slot gets called every time the QNetworkReply has new data.
    // We read all of its new data and write it into the file.
    // That way we use less RAM than when reading it at the finished()
    // signal of the QNetworkReply
    if (file)
        file->write(reply->readAll());
}

void fAppDownloadDialog::updateDownloadProgress(qint64 bytesRead, qint64 totalBytes)
{
    if (httpRequestAborted)
        return;

    // downloadProgressDialog->setMaximum(totalBytes);
    // downloadProgressDialog->setValue(bytesRead);
    QString msg = "Downloading " + this->appName;
    emit updateProgress((int) ((bytesRead * 100) / totalBytes), msg.toStdString(), 100);
}

// When download finished or canceled, this will be called
void fAppDownloadDialog::httpDownloadFinished()
{
  // when canceled
    if (httpRequestAborted) {
        if (file) {
            file->close();
            file->remove();
            delete file;
            file = nullptr;
        }
        reply->deleteLater();
        // downloadProgressDialog->hide();

        ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
        ApplicationPreferences::GetInstance()->SerializePreferences();
        // ApplicationPreferences::GetInstance()->DisplayPreferences();

        // emit cancelDownload();

        return;
    }

    // download finished normally
    // downloadProgressDialog->hide();
    file->flush();
    file->close();

    // get redirection url
    QVariant redirectionTarget = reply->attribute(QNetworkRequest::RedirectionTargetAttribute);
    if (reply->error()) {
        file->remove();
        QMessageBox::information(this, tr("HTTP"),
                                tr("Download failed: %1.")
                                .arg(reply->errorString()));
        // ui->downloadButton->setEnabled(true  );
    } else if (!redirectionTarget.isNull()) {
        QUrl newUrl = url.resolved(redirectionTarget.toUrl());
        if (QMessageBox::question(this, tr("HTTP"),
                                tr("Redirect to %1 ?").arg(newUrl.toString()),
                                QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes) {
            url = newUrl;
            reply->deleteLater();
            file->open(QIODevice::WriteOnly);
            file->resize(0);
            startRequest(url);
            return;
        }
    } else {
        QString fileName = QFileInfo(QUrl(qInputLink).path()).fileName();

    }

    reply->deleteLater();
    reply = nullptr;
    delete file;
    file = nullptr;
    manager = nullptr;

    ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("true").toString());
    ApplicationPreferences::GetInstance()->SerializePreferences();
    // ApplicationPreferences::GetInstance()->DisplayPreferences();

    emit doneDownload(fullPath, extractPath);
}

// During the download progress, it can be canceled
void fAppDownloadDialog::cancelDownload()
{
    httpRequestAborted = true;
    reply->abort();

    ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
    ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
    ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
    ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
    ApplicationPreferences::GetInstance()->SerializePreferences();
    // ApplicationPreferences::GetInstance()->DisplayPreferences();

    // emit cancelDownload();

    this->close();
}