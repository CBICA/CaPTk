#ifndef ApplicationDownloadManager_H
#define ApplicationDownloadManager_H

#include <QObject>
#include <QMutex>
#include <QProgressDialog>

#include "fAppDownloadDialog.h"
#include "yaml-cpp/node/node.h"

class ApplicationDownloadManager : public QObject
{
	Q_OBJECT

public:
	//! constructor/desctrucor
	ApplicationDownloadManager() = default;
	~ApplicationDownloadManager() = default;
	
	QProgressDialog *extractProgressDialog;

	//! setters/getters
	void SetName(QString appName);
	QString GetName() const;

	std::string getApplication(QString appName);
	// std::string getApplicationCLI(QString appName);

private:

	Q_DISABLE_COPY(ApplicationDownloadManager)

	fAppDownloadDialog appDownloadDialog;
	QString m_AppName;
	
	YAML::Node m_appDownloadConfigs = YAML::LoadFile(getCaPTkDataDir() + "/links.yaml");

	void appDownload(bool isCLI);

private slots:
	void startUnzip(QString fullPath, QString extractPath);
	void doneUnzip();

signals:
	void updateProgress(int progress, std::string message, int max);
};

#endif 
