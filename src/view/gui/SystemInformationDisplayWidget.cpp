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
}

void SystemInformationDisplayWidget::SetupUi()
{
	this->setWindowTitle("System Information");
	this->setMinimumWidth(550);

	verticalLayout = new QVBoxLayout(this);
	verticalLayout->setObjectName(QStringLiteral("verticalLayout"));

	horizontalLayout = new QHBoxLayout();
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