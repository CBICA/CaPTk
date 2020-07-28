#include "ApplicationDownloadManager.h"

#include <QSettings>
#include <QFile>
#include <QDebug>
#include <CaPTkDefines.h>

//#include "CaPTkUtils.h"
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

std::string ApplicationDownloadManager::getApplication(QString appName, bool isCLI) {
	this->m_AppName = appName;
	this->isCLI = isCLI;
	m_logger.Write("Application requested: " + appName.toStdString());
    std::string scriptToCall = getApplicationDownloadPath(this->m_AppName.toStdString());

    if(!this->IsApplicationAvailable(appName)) {
		ApplicationPreferences::GetInstance()->DeSerializePreferences();
		bool downloadStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadStartedStatus()).toBool();
		bool downloadFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadFinishedStatus()).toBool();
		bool extractionStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionStartedStatus()).toBool();
		bool extractionFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionFinishedStatus()).toBool();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();

		m_logger.Write(appName.toStdString() + " was not found. Checking for ongoing installation...");

		if(downloadStarted && !downloadFinished)
		{
			if (!this->isCLI) {
				QMessageBox::information(&appDownloadDialog,tr("Install"),"Download is in progress");
			}
			m_logger.Write("An ongoing download of " + appName.toStdString() + " was detected. ");
			qDebug() << "Download is in progress\n";
			emit extractResult(false);

			return "";
		}

		if(extractionStarted && !extractionFinished)
		{
			if (!this->isCLI) {
				QMessageBox::information(&appDownloadDialog ,tr("Extract"),"Extraction is in progress");
			}
			m_logger.Write("An ongoing extraction of " + appName.toStdString() + " was detected. ");
			qDebug() << "Extraction is in progress\n";
			emit extractResult(false);

			return "";
		}
		m_logger.Write("No ongoing installation was detected.");
		m_logger.Write("Preparing to install " + appName.toStdString());
		appDownload();
		
		return "";
	}
	m_logger.Write(appName.toStdString() + " application was installed at: " + scriptToCall );
    return scriptToCall;
}

bool ApplicationDownloadManager::IsApplicationAvailable(QString app)
{
    bool available = false;
    ApplicationPreferences::GetInstance()->DeSerializePreferences();
    bool downloadStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadStartedStatus()).toBool();
    bool downloadFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadFinishedStatus()).toBool();
    bool extractionStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionStartedStatus()).toBool();
    bool extractionFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionFinishedStatus()).toBool();

	std::string scriptToCall = getApplicationDownloadPath(this->m_AppName.toStdString());

    if(downloadStarted && downloadFinished && extractionStarted && extractionFinished && !scriptToCall.empty())
        available = true;
    else {
        available = false;
    }
    return available;
}

void ApplicationDownloadManager::appDownload()
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
	m_logger.Write("Preparing to download " + m_AppName.toStdString() + " from " + downloadLink +
		" into " + downloadFolder);
    appDownloadDialog.SetPaths(downloadFolder);
    appDownloadDialog.SetDownloadLink(downloadLink);
	if (this->isCLI) {
		qDebug() << "Downloading " << this->m_AppName << "\n";

		appDownloadDialog.initDownload();
	}
	else {
		appDownloadDialog.exec();
	}

	connect( &appDownloadDialog, SIGNAL(updateProgress(int, std::string, int)), this, SLOT(updateProgressSlot(int, std::string, int)));   
	connect( &appDownloadDialog, SIGNAL(doneDownload(QString, QString)), this, SLOT(startUnzip(QString, QString)));   
}

void ApplicationDownloadManager::updateProgressSlot(int progress, std::string message, int max) {
	
	// qDebug() << QString::fromStdString(message) << endl;

	emit updateProgressSignal(progress, message, max);
}

void ApplicationDownloadManager::startUnzip(QString fullPath, QString extractPath) 
{
	m_logger.Write("Looking for the downloaded file...");

    if (cbica::isFile(fullPath.toStdString())) {
		ThreadedExtraction* asyncExtract = new ThreadedExtraction();

		connect(asyncExtract, SIGNAL(updateProgressSignal(int, std::string, int)), this, SLOT(updateProgressSlot(int, std::string, int))); 
        connect(asyncExtract, SIGNAL(resultReady(bool)), this, SLOT(doneUnzip(bool)));
		connect(asyncExtract, &ThreadedExtraction::finished, asyncExtract, &QObject::deleteLater);

		asyncExtract->setFullPath(fullPath);
		asyncExtract->setExtractPath(extractPath);
		asyncExtract->setAppName(this->m_AppName);

		qDebug() << "\nExtracting " << this->m_AppName << "\n";
		m_logger.Write("Starting extraction of " + m_AppName.toStdString() + " from " + fullPath.toStdString()
			+ " into " + extractPath.toStdString());

		updateProgressSlot(0, "Extracting " + this->m_AppName.toStdString(), 100);

		asyncExtract->start();
	}
	else
	{
		m_logger.WriteError("Failed to find the archive file at " + fullPath.toStdString()
			+ ". Installation cannot proceed.");
	}
}

void ApplicationDownloadManager::doneUnzip(bool extracted) {
    if(!extracted){
		ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		m_logger.WriteError("Extraction failed. Installation cannot continue.");
		if (!this->isCLI) {
			QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction failed");
		}
		updateProgressSlot(0, "Install " + this->m_AppName.toStdString() + " not completed", 100);

        qDebug() << "Install " << this->m_AppName << " failed. Please retry.\n";

	}
	else {
		m_logger.Write("Extraction finished successfully.");
		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		if (!this->isCLI) {
			QMessageBox::information(&appDownloadDialog, tr("Extraction"),"Extraction done");
		}
		updateProgressSlot(100, "Install " + this->m_AppName.toStdString() + " completed", 100);
		m_logger.Write("Succesfully installed " + m_AppName.toStdString());
		qDebug() << "Install " << this->m_AppName << " completed\n";
		qDebug() << "Please rerun " << this->m_AppName << "\n";	
	}

    emit extractResult(extracted);
}
