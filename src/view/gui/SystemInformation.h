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

	//! Get Basic OS information
	void GetBasicOSInformation();

	//! Get Detailed OS information
	void GetDetailedOSInformation();

	//! Get memory information
	void GetMemoryInformation();

	//! Get OpenGL information
	void GetOpenGLInformation();

	QStringList m_InfoList;

};

#endif // SYSTEMINFORMATION_H
