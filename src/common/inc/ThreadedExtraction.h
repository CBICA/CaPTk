#ifndef ASYNCEXTRACT_H
#define ASYNCEXTRACT_H

#include <QObject>
#include <QThread>
#include <QFile>
#include <QDebug>

#include "QZipReader.h"
#include "ApplicationPreferences.h"

class ThreadedExtraction : public QThread
{
	Q_OBJECT
	void run() override;

public:
	ThreadedExtraction() = default;
	~ThreadedExtraction() = default;

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
	Q_DISABLE_COPY(ThreadedExtraction)

	QString fullPath;
	QString extractPath;
	QString appName;

private slots:
	void updateProgressSlot(int progress);

signals:
	void updateProgressSignal(int progress, std::string message, int max);
    void resultReady(bool);
};

#endif
