#ifndef CAPTKDOCKWIDGET_H
#define CAPTKDOCKWIDGET_H

#include <QDockWidget>


/**
\class CaPTkDockWidget
\brief This customized QDockWidget allows you to propagate drag-n-drop events to other widgets through signals.

*/
class CaPTkDockWidget : public QDockWidget
{
    Q_OBJECT
public:
    explicit CaPTkDockWidget(QWidget *parent = nullptr);

protected:
    virtual void dragEnterEvent(QDragEnterEvent* event) override;
    virtual void dropEvent(QDropEvent* event) override;

	//handler for close
	virtual void closeEvent(QCloseEvent *event) override;

protected slots:
	//! Dock/undock behaviour changed
	void toolTabDockChanged(bool bUnDocked);

signals:
    void dragEnteredDockWidget(QDragEnterEvent*);
    void droppedOnDockWidget(QDropEvent*);

	//signal to close the application
	void close(); 

};

#endif // CAPTKDOCKWIDGET_H
