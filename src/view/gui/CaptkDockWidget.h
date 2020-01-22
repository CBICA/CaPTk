#ifndef CAPTKDOCKWIDGET_H
#define CAPTKDOCKWIDGET_H

#include <QDockWidget>


/**
\class CaptkDockWidget
\brief This customized QDockWidget allows you to propagate drag-n-drop events to other widgets through signals.

*/
class CaptkDockWidget : public QDockWidget
{
    Q_OBJECT
public:
    explicit CaptkDockWidget(QWidget *parent = nullptr);
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);

signals:
    void dragEnteredDockWidget(QDragEnterEvent*);
    void droppedOnDockWidget(QDropEvent*);

};

#endif // CAPTKDOCKWIDGET_H
