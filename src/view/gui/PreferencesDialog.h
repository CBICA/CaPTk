#ifndef PREFERENCESDIALOG_H
#define PREFERENCESDIALOG_H

#include <QDialog>
class QFontDialog;
class QColorDialog;
class QVBoxLayout;
class QWidget;
class QHBoxLayout;
class QListWidget;
class QStackedWidget;
class QDialogButtonBox;

class PreferencesDialog : public QDialog
{
    Q_OBJECT

public:
    explicit PreferencesDialog(QWidget *parent = nullptr);
    ~PreferencesDialog();

    QFontDialog* GetFontDialog();

    void SetApplicationFont();

    QString GetColorStyleSheet();

public slots:

    void OnItemSelectionChanged();

private:

	void SetupUi();
	void retranslateUi(QDialog *PreferencesDialog);

    //Ui::PreferencesDialog *ui;
    QFontDialog *m_FontDialog;
    QColorDialog *m_ColorDialog;
    QString m_colorStyleSheetString;

	QVBoxLayout *verticalLayout;
	QWidget *widget;
	QHBoxLayout *horizontalLayout;
	QListWidget *listWidget;
	QStackedWidget *stackedWidget;
	QDialogButtonBox *buttonBox;
};

#endif // PREFERENCESDIALOG_H
