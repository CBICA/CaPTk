#include "ApplicationDownloadManager.h"

#include <QSettings>
#include <QFile>
#include <QDebug>
#include <CaPTkDefines.h>

#include "cbicaLogging.h"
#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "ApplicationPreferences.h"
#include "ThreadedExtraction.h"

void ApplicationDownloadManager::SetName(QString appName)
{
	this->m_AppName = appName;
}

QString ApplicationDownloadManager::GetName() const
{
	return m_AppName;
}

std::string ApplicationDownloadManager::getApplication(QString appName) {
	this->m_AppName = appName;

	std::string scriptToCall = getApplicationDownloadPath(this->m_AppName.toStdString());

	if (scriptToCall.empty()) {
		ApplicationPreferences::GetInstance()->DeSerializePreferences();
		bool downloadStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadStartedStatus()).toBool();
		bool downloadFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadFinishedStatus()).toBool();
		ApplicationPreferences::GetInstance()->DisplayPreferences();

		if(downloadStarted && !downloadFinished)
		{
			QMessageBox::information(&appDownloadDialog,tr("Download"),"Download in progress");
			return "";
		}

		bool isCLI = false;
		appDownload(isCLI);
		
		return "";
	}
	else {
		ApplicationPreferences::GetInstance()->DeSerializePreferences();
		bool extractionStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionStartedStatus()).toBool();
		bool extractionFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionFinishedStatus()).toBool();
		ApplicationPreferences::GetInstance()->DisplayPreferences();

		if(extractionStarted && !extractionFinished)
		{
			QMessageBox::information(&appDownloadDialog ,tr("Extract"),"Extraction in progress");

			return "";
		}
	}

	return scriptToCall;
}

void ApplicationDownloadManager::appDownload(bool isCLI)
{
	std::string linkyml = "";

	#ifdef _WIN32
	linkyml = "Windows";
	#elif __APPLE__
	linkyml = "macOS";
	#else
	linkyml = "Linux";
	#endif

	std::string downloadLink = m_appDownloadConfigs["appsDownload"][this->m_AppName.toStdString()][linkyml].as<std::string>();
	
	appDownloadDialog.SetPaths(downloadFolder);
	appDownloadDialog.SetDownloadLink(downloadLink);
	if (isCLI) {
		appDownloadDialog.initDownload();
	}
	else {
		appDownloadDialog.exec();
	}

	connect( &appDownloadDialog, SIGNAL(updateProgressDownload(int, std::string, int)), this, SLOT(updateProgressDownload(int, std::string, int)));   
	connect( &appDownloadDialog, SIGNAL(doneDownload(QString, QString)), this, SLOT(startUnzip(QString, QString)));   
}

void ApplicationDownloadManager::updateProgress(int progress, std::string message, int max) {
	emit updateProgressDownload(progress, "Downloading " + this->m_AppName.toStdString(), max);
}


void ApplicationDownloadManager::startUnzip(QString fullPath, QString extractPath) 
{
	if (cbica::isFile(fullPath.toStdString())) {

		ThreadedExtraction* asyncExtract = new ThreadedExtraction();

		connect(asyncExtract, SIGNAL(resultReady(QString)), this, SLOT(doneUnzip()));
		connect(asyncExtract, &ThreadedExtraction::finished, asyncExtract, &QObject::deleteLater);

		asyncExtract->setFullPath(fullPath);
		asyncExtract->setExtractPath(extractPath);
		asyncExtract->setAppName(this->m_AppName);

		asyncExtract->start();

		QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction has started in the background and will take 10-15 minutes");

		// extractProgressDialog = new QProgressDialog(&appDownloadDialog);
    	// extractProgressDialog->setObjectName(QString::fromUtf8("ProgressDialog"));
		// extractProgressDialog->setLabelText(tr("Extracting %1").arg(this->m_AppName));
		// extractProgressDialog->setRange(0, 0);
		// extractProgressDialog->setCancelButton(0);
		emit updateProgressExtract(50, "Installing " + this->m_AppName.toStdString(), 100);
	}
}

void ApplicationDownloadManager::doneUnzip() {

	if (getApplicationDownloadPath(this->m_AppName.toStdString()).empty()) {
		
		// extractProgressDialog->cancel();
		QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction failed");
		// qDebug() << "Extraction failed" << endl;
		ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		emit updateProgressExtract(0, "Install " + this->m_AppName.toStdString() + " not completed", 100);

	}
	else {
		// extractProgressDialog->cancel();
		QMessageBox::information(&appDownloadDialog, tr("Extraction"),"Extraction done");
		// qDebug() << "Extraction done" << endl;

		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		emit updateProgressExtract(100, "Install " + this->m_AppName.toStdString() + " completed", 100);

	}
}