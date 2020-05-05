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
	virtual void closeEvent(QCloseEvent *event) override;

signals:
    void dragEnteredDockWidget(QDragEnterEvent*);
    void droppedOnDockWidget(QDropEvent*);
	void close(); 

};

#endif // CAPTKDOCKWIDGET_H
