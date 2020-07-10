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

	//! Get OS information
	void GetOSInformation();

	QStringList m_InfoList;

};

#endif // SYSTEMINFORMATION_H
