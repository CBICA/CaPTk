#include "StandaloneApp.h"

#include <QSettings>
#include <QFile>
#include <QDebug>
#include <CaPTkDefines.h>

#include "cbicaLogging.h"
#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "AsyncExtract.h"
#include "ApplicationPreferences.h"

// StandaloneApp* StandaloneApp::m_Instance = nullptr;
// QMutex StandaloneApp::m_Mutex;

// StandaloneApp* StandaloneApp::GetInstance()
// {
// 	if (m_Instance == nullptr)
// 	{
// 		m_Mutex.lock();
// 		m_Instance = new StandaloneApp();
// 		// m_Instance.SetName(appName);
// 		m_Mutex.unlock();
// 	}
// 	return m_Instance;
// }

void StandaloneApp::SetName(QString appName)
{
	this->m_AppName = appName;
}

QString StandaloneApp::GetName() const
{
	return m_AppName;
}

std::string StandaloneApp::getStandaloneApp(QString appName) {
	// this = pMainWindow;
	this->m_AppName = appName;

	std::string scriptToCall = getApplicationDownloadPath(this->m_AppName.toStdString());

	// StandaloneApp* stlapps = StandaloneApp::GetInstance();

	// stlapps->RetreiveAppSetting(QString::fromStdString(appName));
	// stlapps->Debug("Function call");

	// if (!(stlapps->GetAction() == "Download" && stlapps->GetStatus() == "Start")) { // if download is not started
	// 	if (scriptToCall.empty()) { // app not found or delete after extraction
	// 		stlapps->StoreAppSetting("", "", QString::fromStdString(appName));
	// 	}

	// 	if (stlapps->GetAction() == "Extract" && stlapps->GetStatus() == "Done") { // if extraction finished
	// 		scriptToCall = getApplicationDownloadPath(appName);

	// 		stlapps->RetreiveAppSetting(QString::fromStdString(appName));
	// 		stlapps->Debug("Path Set");

	// 		return scriptToCall;
	// 	}
	// 	else if (stlapps->GetAction() == "Extract" && stlapps->GetStatus() == "Start") { // if extraction finished
	// 		ShowErrorMessage("The application is being installed");
	// 		updateProgress(50, "Extracting " + appName);
	// 		return "";
	// 	}
	// 	else if (!(stlapps->GetAction() == "Download" && stlapps->GetStatus() == "Done")) { // if download is never started or not done before
	if (scriptToCall.empty()) {
		appDownload();
		
		return "";
	}

	return scriptToCall;
	// 	} 
	// } 
	// else { // download already started
	// 	ShowErrorMessage("The application is being downloaded");
	// 	return "";
	// }
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

	// ShowErrorMessage(downloadLink);

	appDownloadDialog.SetPaths(downloadFolder);
	appDownloadDialog.SetDownloadLink(downloadLink);
	appDownloadDialog.exec();

	connect( &appDownloadDialog, SIGNAL(doneDownload(QString, QString)), this, SLOT(startUnzip(QString, QString))); 
	// connect( &appDownloadDialog, SIGNAL(startDownload()), this, SLOT(startDownload()));    
	// connect( &appDownloadDialog, SIGNAL(cancelDownload()), this, SLOT(cancelDownload()));    
}

void StandaloneApp::startDownload() 
{
	// StandaloneApp* stlapps = StandaloneApp::GetInstance();

	// stlapps->RetreiveAppSetting(appName);
	// stlapps->Debug("Download Start");

	// stlapps->StoreAppSetting("Download", "Start", appName);
}

void StandaloneApp::cancelDownload() 
{
	// StandaloneApp* stlapps = StandaloneApp::GetInstance();

	// stlapps->RetreiveAppSetting(appName);
	// stlapps->Debug("Cancel Download");

	// stlapps->StoreAppSetting("", "", appName);
}

void StandaloneApp::startUnzip(QString fullPath, QString extractPath) 
{
	if (cbica::isFile(fullPath.toStdString())) {
		//  StandaloneApp* stlapps = StandaloneApp::GetInstance();

		//  stlapps->RetreiveAppSetting(appName);
		//  stlapps->Debug("Done download");

		//  stlapps->StoreAppSetting("Download", "Done", appName);

		//  updateProgress(50, "Extracting " + appName.toStdString());
		ASyncExtract* asyncExtract = new ASyncExtract();

		connect(asyncExtract, SIGNAL(resultReady(QString)), this, SLOT(doneUnzip(QString)));
		connect(asyncExtract, &ASyncExtract::finished, asyncExtract, &QObject::deleteLater);

		asyncExtract->setFullPath(fullPath);
		asyncExtract->setExtractPath(extractPath);
		asyncExtract->setAppName(this->m_AppName);

		asyncExtract->start();
	}
}

void StandaloneApp::doneUnzip() {
	//   StandaloneApp* stlapps = StandaloneApp::GetInstance();

	if (getApplicationDownloadPath(this->m_AppName.toStdString()).empty()) {

		// updateProgress(0, "Extracting " + this->m_AppName.toStdString() + " failed");

		// ShowErrorMessage("Installation failed. Please re-run installtion.");
		//  stlapps->RetreiveAppSetting(appName);
		//  stlapps->Debug("Extraction failed");

		//  stlapps->StoreAppSetting("", "", appName);
		// QMessageBox::information(this,tr("Extraction"),"Extraction failed");
		qDebug() << "Extraction failed" << endl;

	}
	else {
		// updateProgress(100, "Extracting " + this->m_AppName.toStdString() + " done");
		// QMessageBox::information(this, tr("Extraction"),"Extraction done");
		qDebug() << "Extraction done" << endl;


		//  stlapps->RetreiveAppSetting(appName);
		//  stlapps->Debug("Extraction done");
		//  stlapps->StoreAppSetting("Extract", "Done", appName);
	}
}