#include "ThreadedDownload.h"
#include <QDebug>
#include "ApplicationDownloadManager.h"

ThreadedDownload::ThreadedDownload(QObject *parent) :
    QThread(parent)
{
}

// We overrides the QThread's run() method here
// run() will be called when a thread starts
// the code will be shared by all threads

void ThreadedDownload::run()
{
    ApplicationDownloadManager* appDownloadMngr = new ApplicationDownloadManager();
    std::string libraPath = appDownloadMngr->getApplication("libra", true);
}