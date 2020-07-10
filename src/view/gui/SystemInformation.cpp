#include "systeminformation.h"
#include <QOperatingSystemVersion>
#include <QStorageInfo>
#include <QNetworkInterface>

SystemInformation::SystemInformation()
{
	this->GetOSInformation();
}

QStringList SystemInformation::GetSystemInformation()
{
	return m_InfoList;
}

void SystemInformation::GetOSInformation()
{
	this->m_InfoList << "****OS Information****";
	this->m_InfoList << "OS Version: " + QOperatingSystemVersion::current().name();
	this->m_InfoList << "Major OS Version: " + QString::number(QOperatingSystemVersion::current().majorVersion());
	this->m_InfoList << "Minor Version: " + QString::number(QOperatingSystemVersion::current().minorVersion());
	this->m_InfoList << "Micro Version: " + QString::number(QOperatingSystemVersion::current().microVersion());

}
