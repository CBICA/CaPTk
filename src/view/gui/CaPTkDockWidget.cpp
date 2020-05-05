#include "CaPTkDockWidget.h"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QMimeData>
#include <iostream>

CaPTkDockWidget::CaPTkDockWidget(QWidget *parent) : QDockWidget(parent)
{
    // We must specifically allow drops on the DockWidget.
    this->setAcceptDrops(true);
}

void CaPTkDockWidget::dragEnterEvent(QDragEnterEvent* event) 
{
    /* This emits the drag-enter event to whatever this dockwidget is connected to.
     * calling event->acceptProposedAction() is the responsibility of the receiver,
     * and this must be performed before any QDropEvent will function.
    **/
    emit dragEnteredDockWidget(event);
}

void CaPTkDockWidget::dropEvent(QDropEvent *event) 
{
    // This emits the accepted drop event. Handling the contents is the responsibility of the receiver.
    emit droppedOnDockWidget(event);
}

void CaPTkDockWidget::closeEvent(QCloseEvent * event)
{
	//we reach here when Alt+F4 is pressed to close the dockwidget

	//we handle this case only when the dockwidget is floating
	if (this->isFloating())
	{
		//we don't want to handle the dockwidget close
		event->ignore();

		//instead we want to close the application 
		//we throw a signal to the mainwindow to close
		emit close();
	}
}
