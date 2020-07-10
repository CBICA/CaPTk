#include "systeminformation.h"
#include <QOperatingSystemVersion>
#include <QStorageInfo>
#include <QNetworkInterface>

SystemInformation::SystemInformation()
{

}

QStringList SystemInformation::GetSystemInformation()
{
	return m_InfoList;
}
