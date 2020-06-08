#include "StandaloneApp.h"

#include <QSettings>
#include <QFile>
#include <QDebug>
#include <CaPTkDefines.h>

#include "cbicaLogging.h"
#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "ApplicationPreferences.h"

void StandaloneApp::SetName(QString appName)
{
	this->m_AppName = appName;
}

QString StandaloneApp::GetName() const
{
	return m_AppName;
}

std::string StandaloneApp::getStandaloneApp(QString appName) {
	this->m_AppName = appName;

	std::string scriptToCall = getApplicationDownloadPath(this->m_AppName.toStdString());

	if (scriptToCall.empty()) {
		appDownload();
		
		return "";
	}
	else {
		ApplicationPreferences::GetInstance()->DeSerializePreferences();
		bool extractionStarted = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionStartedStatus()).toBool();
		bool extractionFinished = QVariant(ApplicationPreferences::GetInstance()->GetLibraExtractionFinishedStatus()).toBool();

		ApplicationPreferences::GetInstance()->DisplayPreferences();

		if(extractionStarted && !extractionFinished)
		{
			QMessageBox::information(this,tr("Extract"),"Extraction in progress");
			this->close();

			return "";
		}
	}

	return scriptToCall;
}

void StandaloneApp::appDownload()
{
	std::string linkyml = "";

	#ifdef _WIN32
	linkyml = "Windows";
	#elif __APPLE__
	linkyml = "macOS";
	#else
	linkyml = "Linux";
	#endif

	std::string downloadLink = m_appDownloadConfigs["apps"][this->m_AppName.toStdString()][linkyml].as<std::string>();

	appDownloadDialog.SetPaths(downloadFolder);
	appDownloadDialog.SetDownloadLink(downloadLink);
	appDownloadDialog.exec();

	connect( &appDownloadDialog, SIGNAL(doneDownload(QString, QString)), this, SLOT(startUnzip(QString, QString)));   
}

void StandaloneApp::startUnzip(QString fullPath, QString extractPath) 
{
	if (cbica::isFile(fullPath.toStdString())) {

		ASyncExtract* asyncExtract = new ASyncExtract();

		connect(asyncExtract, SIGNAL(resultReady(QString)), this, SLOT(doneUnzip()));
		connect(asyncExtract, &ASyncExtract::finished, asyncExtract, &QObject::deleteLater);

		asyncExtract->setFullPath(fullPath);
		asyncExtract->setExtractPath(extractPath);
		asyncExtract->setAppName(this->m_AppName);

		asyncExtract->start();
	}
}

void StandaloneApp::doneUnzip() {

	if (getApplicationDownloadPath(this->m_AppName.toStdString()).empty()) {

		QMessageBox::information(NULL,tr("Extraction"),"Extraction failed");
		// qDebug() << "Extraction failed" << endl;
		ApplicationPreferences::GetInstance()->SetLibraDownloadStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraDownloadFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();
	}
	else {
		QMessageBox::information(NULL, tr("Extraction"),"Extraction done");
		// qDebug() << "Extraction done" << endl;

		ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();
	}
}