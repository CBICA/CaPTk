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

	connect( &appDownloadDialog, SIGNAL(updateProgress(int, std::string, int)), this, SLOT(updateProgressSlot(int, std::string, int)));   
	connect( &appDownloadDialog, SIGNAL(doneDownload(QString, QString)), this, SLOT(startUnzip(QString, QString)));   
}

void ApplicationDownloadManager::updateProgressSlot(int progress, std::string message, int max) {
	if (message.find("Download") == std::string::npos) {
		QMessageBox::information(&appDownloadDialog, tr("Info"), QString::fromStdString(message));
	}

	qDebug() << QString::fromStdString(message) << endl;

	emit updateProgressSignal(progress, message, max);
}

void ApplicationDownloadManager::startUnzip(QString fullPath, QString extractPath) 
{
	if (cbica::isFile(fullPath.toStdString())) {

		ThreadedExtraction* asyncExtract = new ThreadedExtraction();

		connect(asyncExtract, SIGNAL(updateProgressSignal(int, std::string, int)), this, SLOT(updateProgressSlot(int, std::string, int))); 
		connect(asyncExtract, SIGNAL(resultReady(QString)), this, SLOT(doneUnzip()));
		connect(asyncExtract, &ThreadedExtraction::finished, asyncExtract, &QObject::deleteLater);

		asyncExtract->setFullPath(fullPath);
		asyncExtract->setExtractPath(extractPath);
		asyncExtract->setAppName(this->m_AppName);

		// updateProgressSlot(0, "Installing " + this->m_AppName.toStdString(), 100);

		asyncExtract->start();

		// QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction has started in the background and will take 10-15 minutes");
	}
}

void ApplicationDownloadManager::doneUnzip() {

	if (getApplicationDownloadPath(this->m_AppName.toStdString()).empty()) {
		
		// qDebug() << "Extraction failed" << endl;
		ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		// updateProgressSlot(0, "Install " + this->m_AppName.toStdString() + " not completed", 100);
		QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction failed");
	}
	else {
		// qDebug() << "Extraction done" << endl;

		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		// updateProgressSlot(100, "Install " + this->m_AppName.toStdString() + " completed", 100);
		QMessageBox::information(&appDownloadDialog, tr("Extraction"),"Extraction done");
	}
}