#ifndef ApplicationDownloadManager_H
#define ApplicationDownloadManager_H

#include <QObject>
#include <QMutex>
#include <QProgressDialog>

#include "CaPTkGUIUtils.h"
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

	std::string getApplication(QString appName, bool isCLI);
	// std::string getApplicationCLI(QString appName);

    // void SetDownloadFolder(QString folder)
    // {
    //     this->downloadFolder = folder;
    // }

    // void SetDownloadLink(QString link)
    // {
    //     this->downloadLink = link;
    // }

    bool IsApplicationAvailable(QString app);

private:

	Q_DISABLE_COPY(ApplicationDownloadManager)

	fAppDownloadDialog appDownloadDialog;
	QString m_AppName;
	bool isCLI;
    // QString downloadFolder, downloadLink;
	
    YAML::Node m_appDownloadConfigs = YAML::LoadFile(getCaPTkDataDir() + "/links.yaml");
	cbica::Logging m_logger = cbica::Logging(loggerFile, "");

	void appDownload();

private slots:
	void updateProgressSlot(int progress, std::string message, int max);
	void startUnzip(QString fullPath, QString extractPath);
    void doneUnzip(bool);

signals:
	void updateProgressSignal(int progress, std::string message, int max);
    void extractResult(bool extracted);
};

#endif 
