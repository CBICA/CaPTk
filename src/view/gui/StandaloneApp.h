#ifndef StandaloneApp_H
#define StandaloneApp_H

#include <QObject>
#include <QMutex>
#include <QThread>

#include "fAppDownloadDialog.h"
#include "fMainWindow.h"
#include "yaml-cpp/node/node.h"

class StandaloneApp : public QObject
{
public:
	//! constructor/desctrucor
	StandaloneApp() = default;
	~StandaloneApp() = default;

	//! setters/getters
	void SetName(QString appName);
	QString GetName() const;

	std::string getStandaloneApp(QString appName, fMainWindow* pMainWindow);

private:

	Q_DISABLE_COPY(StandaloneApp)

	//! ivars
	// static StandaloneApp* m_Instance;
	// static QMutex m_Mutex;

	static fMainWindow* m_fMainWindow;
	fAppDownloadDialog appDownloadDialog;
	QString m_AppName;
	
	YAML::Node m_appDownloadConfigs = YAML::LoadFile(getCaPTkDataDir() + "/appsDownloadConfigs.yaml");

	void appDownload();
	void startDownload();
	void cancelDownload();
	void startUnzip(QString fullPath, QString extractPath);
	void doneUnzip();
};

#endif 
