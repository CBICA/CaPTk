#include "systeminformation.h"
#include <QOperatingSystemVersion>
#include <QStorageInfo>
#include <QNetworkInterface>

#include "CaPTkUtils.h"

SystemInformation::SystemInformation()
{
	this->GetBasicOSInformation();
	this->GetDetailedOSInformation();
	this->GetMemoryInformation();
}

QStringList SystemInformation::GetSystemInformation()
{
	return m_InfoList;
}

void SystemInformation::GetBasicOSInformation()
{
	this->m_InfoList << "****Basic OS Information****";
	this->m_InfoList << "OS Version: " + QOperatingSystemVersion::current().name();
	this->m_InfoList << "Major OS Version: " + QString::number(QOperatingSystemVersion::current().majorVersion());
	this->m_InfoList << "Minor Version: " + QString::number(QOperatingSystemVersion::current().minorVersion());
	this->m_InfoList << "Micro Version: " + QString::number(QOperatingSystemVersion::current().microVersion());
	this->m_InfoList << "\n";
}

void SystemInformation::GetDetailedOSInformation()
{
	QSysInfo systemInfo;
	this->m_InfoList << "****Detailed OS Information****";
	this->m_InfoList << "Cpu Architecture: " + systemInfo.currentCpuArchitecture();
	this->m_InfoList << "Kernel Type: " + systemInfo.kernelType();
	this->m_InfoList << "Kernel Version: " + systemInfo.kernelVersion();
	this->m_InfoList << "Machine Host Name: " + systemInfo.machineHostName();
	this->m_InfoList << "Product Type: " + systemInfo.productType();
	this->m_InfoList << "Product Version: " + systemInfo.productVersion();
	this->m_InfoList << "Byte Order: " + systemInfo.buildAbi();
	this->m_InfoList << "Pretty ProductName: " + systemInfo.prettyProductName();
	this->m_InfoList << "\n";
}

void SystemInformation::GetMemoryInformation()
{
	/**** Get total amount of ram ****/
	unsigned long long ram = getTotalInstalledMemory();
	this->m_InfoList << "****Memory Information****";
	this->m_InfoList << "RAM: " + QString::number(ram / (1024.0 * 1024 * 1024), 'f', 2) + "GB";
	this->m_InfoList << "\n";
}
