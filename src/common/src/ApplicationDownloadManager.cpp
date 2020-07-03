#include "ApplicationDownloadManager.h"

#include <QSettings>
#include <QFile>
#include <QDebug>
#include <CaPTkDefines.h>

//#include "cbicaLogging.h"
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

    std::string scriptToCall = getApplicationDownloadPath(this->m_AppName.toStdString());

    // std::string scriptToCall = this->downloadFolder.toStdString();

    if(!this->IsApplicationAvailable(appName)){
		ApplicationPreferences::GetInstance()->DeSerializePreferences();
		bool downloadStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadStartedStatus()).toBool();
		bool downloadFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraDownloadFinishedStatus()).toBool();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();

		if(downloadStarted && !downloadFinished)
		{
			if (!this->isCLI) {
				QMessageBox::information(&appDownloadDialog,tr("Download"),"Download in progress");
			}
			std::cout << "Download in progress\n";

			return "";
		}

		appDownload();
		
		return "";
	}
	else {
		ApplicationPreferences::GetInstance()->DeSerializePreferences();
		bool extractionStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionStartedStatus()).toBool();
		bool extractionFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionFinishedStatus()).toBool();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();

		if(extractionStarted && !extractionFinished)
		{
			if (!this->isCLI) {
				QMessageBox::information(&appDownloadDialog ,tr("Extract"),"Extraction in progress");
			}
			std::cout << "Extraction in progress\n";

			return "";
		}
	}

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

    if(downloadStarted && downloadFinished && extractionStarted && extractionFinished)
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
	
    appDownloadDialog.SetPaths(downloadFolder);
    appDownloadDialog.SetDownloadLink(downloadLink);
	if (this->isCLI) {
		std::cout << "Downloading " << this->m_AppName.toStdString() << ": ";

		appDownloadDialog.initDownload();
	}
	else {
		appDownloadDialog.exec();
	}

	connect( &appDownloadDialog, SIGNAL(updateProgress(int, std::string, int)), this, SLOT(updateProgressSlot(int, std::string, int)));   
	connect( &appDownloadDialog, SIGNAL(doneDownload(QString, QString)), this, SLOT(startUnzip(QString, QString)));   
}

void ApplicationDownloadManager::updateProgressSlot(int progress, std::string message, int max) {
	
//	qDebug() << QString::fromStdString(message) << endl;

	emit updateProgressSignal(progress, message, max);
}

void ApplicationDownloadManager::startUnzip(QString fullPath, QString extractPath) 
{
    if (cbica::isFile(fullPath.toStdString())) {

		ThreadedExtraction* asyncExtract = new ThreadedExtraction();

		connect(asyncExtract, SIGNAL(updateProgressSignal(int, std::string, int)), this, SLOT(updateProgressSlot(int, std::string, int))); 
        connect(asyncExtract, SIGNAL(resultReady(bool)), this, SLOT(doneUnzip(bool)));
		connect(asyncExtract, &ThreadedExtraction::finished, asyncExtract, &QObject::deleteLater);

		asyncExtract->setFullPath(fullPath);
		asyncExtract->setExtractPath(extractPath);
		asyncExtract->setAppName(this->m_AppName);

		std::cout << "\nInstalling " << this->m_AppName.toStdString() << ": ";
		// qDebug() << "Extracting" << endl;

		updateProgressSlot(0, "Installing " + this->m_AppName.toStdString(), 100);

		asyncExtract->start();

//        if(this->isCLI)
//            asyncExtract->wait();
		// QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction has started in the background");
	}
}

void ApplicationDownloadManager::doneUnzip(bool extracted) {

    //if (getApplicationDownloadPath(this->m_AppName.toStdString()).empty()) {
		
    if(!extracted){
		// qDebug() << "Extraction failed" << endl;
		ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		updateProgressSlot(0, "Install " + this->m_AppName.toStdString() + " not completed", 100);
		if (!this->isCLI) {
			QMessageBox::information(&appDownloadDialog,tr("Extraction"),"Extraction failed");
		}
        std::cout << "Install " << this->m_AppName.toStdString() << " failed. Please retry.\n";

	}
	else {
		// qDebug() << "Extraction done" << endl;

		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		// ApplicationPreferences::GetInstance()->DisplayPreferences();
		
		updateProgressSlot(100, "Install " + this->m_AppName.toStdString() + " completed", 100);
		if (!this->isCLI) {
			QMessageBox::information(&appDownloadDialog, tr("Extraction"),"Extraction done");
		}
		std::cout << "Install " << this->m_AppName.toStdString() << " completed\n";
		std::cout << "Please rerun " << this->m_AppName.toStdString() << "\n";	
	}

    emit extractResult(extracted);
}
