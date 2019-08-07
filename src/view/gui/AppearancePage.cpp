#include "AppearancePage.h"
#include "ui_AppearancePage.h"

AppearancePage::AppearancePage(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::AppearancePage)
{
    ui->setupUi(this);
}

AppearancePage::~AppearancePage()
{
    delete ui;
}
