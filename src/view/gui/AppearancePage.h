#ifndef APPEARANCEPAGE_H
#define APPEARANCEPAGE_H

#include <QWidget>

namespace Ui {
class AppearancePage;
}

class AppearancePage : public QWidget
{
    Q_OBJECT

public:
    explicit AppearancePage(QWidget *parent = nullptr);
    ~AppearancePage();

private:
    Ui::AppearancePage *ui;
};

#endif // APPEARANCEPAGE_H
