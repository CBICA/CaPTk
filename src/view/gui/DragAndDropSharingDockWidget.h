#ifndef DRAGANDDROPSHARINGDOCKWIDGET_H
#define DRAGANDDROPSHARINGDOCKWIDGET_H

#include <QDockWidget>


/**
\class DragAndDropSharingDockWidget
\brief This customized QDockWidget allows you to propagate drag-n-drop events to other widgets through signals.

*/
class DragAndDropSharingDockWidget : public QDockWidget
{
    Q_OBJECT
public:
    explicit DragAndDropSharingDockWidget(QWidget *parent = nullptr);
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);

signals:
    void dragEnteredDockWidget(QDragEnterEvent*);
    void droppedOnDockWidget(QDropEvent*);

};

#endif // DRAGANDDROPSHARINGDOCKWIDGET_H
