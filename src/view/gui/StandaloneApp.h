#ifndef StandaloneApp_H
#define StandaloneApp_H

#include <QObject>
#include <QMutex>
#include <QThread>

#include "fAppDownloadDialog.h"
#include "yaml-cpp/node/node.h"

class StandaloneApp : public QObject
{
	Q_OBJECT

public:
	//! constructor/desctrucor
	StandaloneApp() = default;
	~StandaloneApp() = default;

	//! setters/getters
	void SetName(QString appName);
	QString GetName() const;

	std::string getStandaloneApp(QString appName);

private:

	Q_DISABLE_COPY(StandaloneApp)

	fAppDownloadDialog appDownloadDialog;
	QString m_AppName;
	
	YAML::Node m_appDownloadConfigs = YAML::LoadFile(getCaPTkDataDir() + "/appsDownloadConfigs.yaml");

	void appDownload();

private slots:
	void startUnzip(QString fullPath, QString extractPath);
	void doneUnzip();
};

#endif 
