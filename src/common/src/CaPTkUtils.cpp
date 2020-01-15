#include "CaPTkUtils.h"

unsigned long long getTotalInstalledMemory()
{
#ifdef _WIN32
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  return status.ullTotalPhys;
#else
  long pages = sysconf(_SC_PHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  return pages * page_size;
#endif
}

bool isSizeOfLoadedFilesTooBig(QStringList files, std::string loggerFile, 
                               float maxPercentage)
{
  // Find total size of all files
  unsigned long long imagesSize = 0;
  for (QString& file : files)
  {
    QFile qFile(file);
    if (qFile.open(QIODevice::ReadOnly)){
        imagesSize += qFile.size();  //when file does open. This size is in bytes
        qFile.close();
    }
  }

  /**** Get total amount of ram ****/
  unsigned long long availableMemory = getTotalInstalledMemory();

  // Log values
  if (loggerFile != "")
  {
    cbica::Logging(loggerFile, 
      "Images size: " + std::to_string(imagesSize));
    cbica::Logging(loggerFile, 
      "Total RAM: " + std::to_string(availableMemory));
  }

  // Compare (maxPercentage is arbitrary, default=0.05 it allows images up to 400MB for a 8GB system)
  if (imagesSize > maxPercentage*availableMemory)
  {
    return true;
  }
  else
  {
    return false;
  }
}