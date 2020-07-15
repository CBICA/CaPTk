#include "SystemInformationDisplayWidget.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QTextEdit>
#include <QPushButton>

SystemInformationDisplayWidget::SystemInformationDisplayWidget(QWidget *parent) :
    QWidget(parent)
{
	this->SetupUi();

	//connections
	connect(this->copyToClipboardPushButton, SIGNAL(clicked()), this, SLOT(OnCopyToClipboardBtnClicked()));
}

SystemInformationDisplayWidget::~SystemInformationDisplayWidget()
{
	if (this->label)
		delete this->label;
	if (this->textEdit)
		delete this->textEdit;
	if (this->copyToClipboardPushButton)
		delete this->copyToClipboardPushButton;
	if (this->horizontalLayout)
		delete this->horizontalLayout;
	if (this->verticalLayout)
		delete this->verticalLayout;
}

void SystemInformationDisplayWidget::SetupUi()
{
	this->setWindowTitle("System Information");
	this->setMinimumWidth(550);

	verticalLayout = new QVBoxLayout(this);
	verticalLayout->setObjectName(QStringLiteral("verticalLayout"));

	horizontalLayout = new QHBoxLayout(this);
	horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));

	copyToClipboardPushButton = new QPushButton("Copy to clipboard", this);
	copyToClipboardPushButton->setObjectName(QStringLiteral("CopyToClipboardPushButton"));

	horizontalLayout->addWidget(copyToClipboardPushButton);

	verticalLayout->addLayout(horizontalLayout);

	textEdit = new QTextEdit(this);
	textEdit->setObjectName(QStringLiteral("textEdit"));
	textEdit->setMinimumHeight(500);

	verticalLayout->addWidget(textEdit);

	label = new QLabel(this);
	label->setObjectName(QStringLiteral("label"));
	label->setText("Please review & edit the above information if needed.");

	verticalLayout->addWidget(label);
}

void SystemInformationDisplayWidget::SetInformation(QStringList sl)
{
	foreach(QString str, sl)
	{
		this->textEdit->append(str);
	}
}

void SystemInformationDisplayWidget::OnCopyToClipboardBtnClicked()
{
	//what to copy?
	this->textEdit->selectAll();

	this->textEdit->copy();

}

void SystemInformationDisplayWidget::ClearInformation()
{
	this->textEdit->clear();
}
