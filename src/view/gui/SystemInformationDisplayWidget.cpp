#include "SystemInformationDisplayWidget.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QTextEdit>
#include <QPushButton>

SystemInformationDisplayWidget::SystemInformationDisplayWidget(QWidget *parent) :
    QWidget(parent)
{
}

void SystemInformationDisplayWidget2::SetupUi()
{
	this->setWindowTitle("System Information");
	this->setMinimumWidth(550);

	verticalLayout = new QVBoxLayout(this);
	verticalLayout->setObjectName(QStringLiteral("verticalLayout"));

	horizontalLayout = new QHBoxLayout();
	horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));

	CopyToClipboardPushButton = new QPushButton("Copy to clipboard", this);
	CopyToClipboardPushButton->setObjectName(QStringLiteral("CopyToClipboardPushButton"));

	horizontalLayout->addWidget(CopyToClipboardPushButton);

	verticalLayout->addLayout(horizontalLayout);

	textEdit = new QTextEdit(this);
	textEdit->setObjectName(QStringLiteral("textEdit"));
	textEdit->setMinimumHeight(500);

	verticalLayout->addWidget(textEdit);

	label = new QLabel(this);
	label->setObjectName(QStringLiteral("label"));
	label->setText("Please edit the above information as you see fit");

	verticalLayout->addWidget(label);

}