#ifndef ApplicationDownloadManager_H
#define ApplicationDownloadManager_H

#include <QObject>
#include <QMutex>

#include "fAppDownloadDialog.h"
#include "yaml-cpp/node/node.h"

class ApplicationDownloadManager : public QObject
{
	Q_OBJECT

public:
	//! constructor/desctrucor
	ApplicationDownloadManager() = default;
	~ApplicationDownloadManager() = default;

	//! setters/getters
	void SetName(QString appName);
	QString GetName() const;

	std::string getApplicationDownloadManager(QString appName);

private:

	Q_DISABLE_COPY(ApplicationDownloadManager)

	fAppDownloadDialog appDownloadDialog;
	QString m_AppName;
	
	YAML::Node m_appDownloadConfigs = YAML::LoadFile(getCaPTkDataDir() + "/links.yaml");

	void appDownload();

private slots:
	void startUnzip(QString fullPath, QString extractPath);
	void doneUnzip();
};

#endif 
