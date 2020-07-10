#ifndef SYSTEMINFORMATION_H
#define SYSTEMINFORMATION_H

#include <QStringList>

class SystemInformation
{
public:
    SystemInformation();
    ~SystemInformation() = default;

	QStringList GetSystemInformation();

private:

	QStringList m_InfoList;

};

#endif // SYSTEMINFORMATION_H
