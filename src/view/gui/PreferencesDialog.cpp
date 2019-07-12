#include "preferencesdialog.h"
//#include "ui_preferencesdialog.h"

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

#include <QFontDialog>
#include <QColorDialog>
#include <QSizePolicy>
#include <QDebug>

PreferencesDialog::PreferencesDialog(QWidget *parent) :
    QDialog(parent)
{
    this->SetupUi();

    //this->m_FontDialog = new QFontDialog(ui->stackedWidget);
    ////this->m_FontDialog->setStyleSheet("background-color:blue;");
    //this->m_FontDialog->setWindowFlags(Qt::SubWindow);
    //this->m_FontDialog->setOptions(
    //            /* do not use native dialog */
    //            QFontDialog::DontUseNativeDialog
    //            /* you don't need to set it, but if you don't set this
    //                                    the "OK" and "Cancel" buttons will show up, I don't
    //                                    think you'd want that. */
    //            | QFontDialog::NoButtons
    //            );
    //this->m_FontDialog->setSizeGripEnabled(false);

    //this->m_ColorDialog = new QColorDialog(ui->stackedWidget);
    ////this->m_ColorDialog->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
    ////this->m_ColorDialog->setStyleSheet("background-color:green;");
    //this->m_ColorDialog->setWindowFlags(Qt::Widget);
    //this->m_ColorDialog->setOptions(
    //            /* do not use native dialog */
    //            QColorDialog::DontUseNativeDialog
    //            /* you don't need to set it, but if you don't set this
    //                                    the "OK" and "Cancel" buttons will show up, I don't
    //                                    think you'd want that. */
    //            | QColorDialog::NoButtons
    //            );
    //ui->stackedWidget->insertWidget(0,this->m_FontDialog);
    //ui->stackedWidget->insertWidget(1,this->m_ColorDialog);
    //ui->stackedWidget->insertWidget(2,new QWidget(ui->stackedWidget));

    //connect(ui->listWidget,SIGNAL(itemSelectionChanged()),this,SLOT(OnItemSelectionChanged()));

}

PreferencesDialog::~PreferencesDialog()
{
    //delete ui;
}

QFontDialog *PreferencesDialog::GetFontDialog()
{
    return this->m_FontDialog;
}

void PreferencesDialog::SetApplicationFont()
{
    //qDebug() << " SetApplicationFont " << endl;
    //qDebug() << ui->stackedWidget->currentIndex() << endl;
    ////if(ui->stackedWidget->currentIndex() == 0)
    //{
    //    qDebug() << " inside if block " << endl;
    //    //qApp->setFont(this->m_FontDialog->currentFont());
    //    QColor c = this->m_ColorDialog->currentColor();
    //    int r = c.red();
    //    int g = c.green();
    //    int b = c.blue();
    //    QString rgbstr = "(" + QString::number(r) + ","
    //            + QString::number(g) + ","
    //            + QString::number(b) + ");";
    //    this->m_colorStyleSheetString = "background-color:rgb" + rgbstr;
    //    qDebug() << " color str = " << m_colorStyleSheetString << endl;
    //    //qApp->setStyleSheet(m_colorStyleSheetString);
    //    //qApp->setProperty("urgent", true);
    //    //ui->listWidget->setStyleSheet(colorstr);
    //}
}

QString PreferencesDialog::GetColorStyleSheet()
{
    return this->m_colorStyleSheetString;
}

void PreferencesDialog::OnItemSelectionChanged()
{
    //int currentItemIndex = ui->listWidget->currentRow();
    //ui->stackedWidget->setCurrentIndex(currentItemIndex);
}

void PreferencesDialog::SetupUi()
{
	if (this->objectName().isEmpty())
		this->setObjectName(QStringLiteral("PreferencesDialog"));
	this->resize(400, 300);
	verticalLayout = new QVBoxLayout(this);
	verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
	widget = new QWidget(this);
	widget->setObjectName(QStringLiteral("widget"));
	widget->setMinimumSize(QSize(0, 0));
	horizontalLayout = new QHBoxLayout(widget);
	horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
	listWidget = new QListWidget(widget);
	new QListWidgetItem(listWidget);

	listWidget->setObjectName(QStringLiteral("listWidget"));
	listWidget->setMinimumSize(QSize(133, 0));
	listWidget->setMaximumSize(QSize(16777215, 16777215));

	horizontalLayout->addWidget(listWidget);

	stackedWidget = new QStackedWidget(widget);
	stackedWidget->setObjectName(QStringLiteral("stackedWidget"));
	stackedWidget->setStyleSheet(QStringLiteral(""));

	horizontalLayout->addWidget(stackedWidget);


	verticalLayout->addWidget(widget);

	buttonBox = new QDialogButtonBox(this);
	buttonBox->setObjectName(QStringLiteral("buttonBox"));
	buttonBox->setOrientation(Qt::Horizontal);
	buttonBox->setStandardButtons(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);

	verticalLayout->addWidget(buttonBox);


	retranslateUi(this);
	QObject::connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	QObject::connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

	stackedWidget->setCurrentIndex(-1);


	QMetaObject::connectSlotsByName(this);
}

void PreferencesDialog::retranslateUi(QDialog *PreferencesDialog)
{
	this->setWindowTitle(QApplication::translate("PreferencesDialog", "Dialog", nullptr));

	const bool __sortingEnabled = listWidget->isSortingEnabled();
	listWidget->setSortingEnabled(false);
	QListWidgetItem *___qlistwidgetitem = listWidget->item(0);
	___qlistwidgetitem->setText(QApplication::translate("PreferencesDialog", "Fonts", nullptr));
	listWidget->setSortingEnabled(__sortingEnabled);
}