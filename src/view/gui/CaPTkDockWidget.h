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
    void dragEnterEvent(QDragEnterEvent* event);
    void dropEvent(QDropEvent* event);
	void closeEvent(QCloseEvent *event);

signals:
    void dragEnteredDockWidget(QDragEnterEvent*);
    void droppedOnDockWidget(QDropEvent*);
	void close(); 

private:
	QWidget* m_parent;

};

#endif // CAPTKDOCKWIDGET_H
