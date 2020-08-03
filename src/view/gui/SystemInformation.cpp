#include "SystemInformation.h"
#include <QOperatingSystemVersion>
#include <QStorageInfo>
#include <QNetworkInterface>

#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOffscreenSurface>

#include "CaPTkUtils.h"

SystemInformation::SystemInformation()
{
	this->GetBasicOSInformation();
	this->GetDetailedOSInformation();
	this->GetMemoryInformation();

	//TBD: The following method needs to be tested on various platforms 
	//If results are incorrect, can be commented till we fix
	this->GetOpenGLInformation(); 
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
	this->m_InfoList << "CPU ABI: " + systemInfo.buildAbi();
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

void SystemInformation::GetOpenGLInformation()
{
	QOffscreenSurface surf;
	surf.create();

	QOpenGLContext ctx;
	ctx.create();
	ctx.makeCurrent(&surf);

	std::string glVersion = std::string(reinterpret_cast<const char*>(ctx.functions()->glGetString(GL_VERSION)));
	std::string glRenderer = std::string(reinterpret_cast<const char*>(ctx.functions()->glGetString(GL_RENDERER)));
	std::string glvendor = std::string(reinterpret_cast<const char*>(ctx.functions()->glGetString(GL_VENDOR)));

	this->m_InfoList << "****OpenGL Information****";
	this->m_InfoList << "OpenGL Vendor: " + QString(glvendor.c_str());
	this->m_InfoList << "OpenGL Version: " + QString(glVersion.c_str());
	this->m_InfoList << "OpenGL Renderer: " + QString(glRenderer.c_str());
	this->m_InfoList << "\n";
}
