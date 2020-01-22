#include "CaptkDockWidget.h"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <iostream>

CaptkDockWidget::CaptkDockWidget(QWidget *parent) : QDockWidget(parent)
{
    // We must specifically allow drops on the DockWidget.
    this->setAcceptDrops(true);
}

void CaptkDockWidget::dragEnterEvent(QDragEnterEvent* event) {
    /* This emits the drag-enter event to whatever this dockwidget is connected to.
     * calling event->acceptProposedAction() is the responsibility of the receiver,
     * and this must be performed before any QDropEvent will function.
    **/
    emit dragEnteredDockWidget(event);
}

void CaptkDockWidget::dropEvent(QDropEvent *event) {
    // This emits the accepted drop event. Handling the contents is the responsibility of the receiver.
    emit droppedOnDockWidget(event);
}