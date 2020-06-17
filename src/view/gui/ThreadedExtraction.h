#ifndef ASYNCEXTRACT_H
#define ASYNCEXTRACT_H

#include <QObject>
#include <QThread>
#include <QFile>
#include <QDebug>

#include "QZipReader.h"
#include "ApplicationPreferences.h"

class ExtractInBackground : public QThread
{
	Q_OBJECT
	void run() override;

public:
	ExtractInBackground() = default;
	~ExtractInBackground() = default;

	void setFullPath(QString fullPath) {
		this->fullPath = fullPath;
	}

	void setExtractPath(QString extractPath) {
		this->extractPath = extractPath;
	}

	void setAppName(QString appName) {
		this->appName = appName;
	}

private:
	Q_DISABLE_COPY(ExtractInBackground)

	QString fullPath;
	QString extractPath;
	QString appName;

signals:
    void resultReady(QString appName);
};

#endif