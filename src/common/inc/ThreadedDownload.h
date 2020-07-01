#ifndef THREADEDDOWNLOAD_H
#define THREADEDDOWNLOAD_H

#include <QThread>

class ThreadedDownload : public QThread
{
    Q_OBJECT
public:
    explicit ThreadedDownload(QObject *parent = 0);
    void run();
signals:
    
public slots:
    
};

#endif // THREADEDDOWNLOAD_H