#ifndef SYSTEMINFORMATIONDISPLAYWIDGET_H
#define SYSTEMINFORMATIONDISPLAYWIDGET_H

#include <QWidget>

class QVBoxLayout;
class QHBoxLayout;
class QPushButton;
class QTextEdit;
class QLabel;

class SystemInformationDisplayWidget : public QWidget
{
    Q_OBJECT
public:
    explicit SystemInformationDisplayWidget(QWidget *parent = nullptr);

	//! set information to display on widget
	void SetInformation(QStringList sl);

public slots:
	void OnCopyToClipboardBtnClicked();

private:

	//set up UI
	void SetupUi();

	//ivars
	QVBoxLayout *verticalLayout;
	QHBoxLayout *horizontalLayout;
	QPushButton *copyToClipboardPushButton;
	QTextEdit *textEdit;
	QLabel *label;
};

#endif // SYSTEMINFORMATIONDISPLAYWIDGET_H
