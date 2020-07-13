#pragma once

#include <cmath>
#include <string>
#include <QMetaType>
#include <QStringList>
#include <QFile>
#include <QIODevice>
#include "CaPTkEnums.h"
#include "CaPTkDefines.h"
#include "cbicaLogging.h"
#include "cbicaUtilities.h"

// For getting the total amount of installed ram
#ifdef _WIN32
#include <windows.h>
#else
// For both linux + mac
#include <unistd.h>
#endif

//! Structure to define a point value to check if it is defined in the image or not
struct PointVal
{
  PointVal()
  {
    x = y = z = value = 0;
  }
  PointVal(int _x, int _y, int _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }
  PointVal(int _x, int _y, int _z, int _value)
  {
    x = _x;
    y = _y;
    z = _z;
    value = _value;
  }
  bool isWithinRange(const int* dims)
  {
    if (x >= 0 && x < dims[0] && y >= 0 && y < dims[1] && z >= 0 && z < dims[2])
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  bool isValid()
  {
    if (x >= 0 && y >= 0 && z >= 0)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  PointVal getInvalidPt()
  {
    return PointVal(-1, -1, -1, -1);
  }
  bool operator==(const PointVal& rhs) const
  {
    if ((x == rhs.x) && (y == rhs.y) && (z == rhs.z))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  float getDistanceFrom(const PointVal& point2)
  {
    return std::sqrt(std::pow(point2.z - z, 2) + std::pow(point2.y - y, 2) + std::pow(point2.x - x, 2));
  }

  int x, y, z, value;
};

Q_DECLARE_METATYPE(PointVal);
Q_DECLARE_METATYPE(std::vector< PointVal >);

struct NonNativeApp
{
  std::string name;
  std::string path;

  //! Default Constructor
  NonNativeApp()
  {
    name = "";
    path = "";
  }

  //! Default Constructor with values
  NonNativeApp(const std::string &inputName, const std::string &inputPath) :
    name(inputName), path(inputPath) {};
};

/**
\brief Guess Image Type

\param str String to guess
\param bool If true, check that the file exists (only succeeds for full paths) and exit if failed. If false, don't check.
\return deduced type
*/
inline int guessImageType(const std::string &fileName, bool checkFileExists = true)
{
	std::string basename = cbica::getFilenameBase(fileName, checkFileExists); 
  int ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;
  std::string fileName_wrap = basename;
  std::transform(fileName_wrap.begin(), fileName_wrap.end(), fileName_wrap.begin(), ::tolower);
  if ((fileName_wrap.find("_t1ce") != std::string::npos) || (fileName_wrap.find("_t1-gad") != std::string::npos) ||
    (fileName_wrap.find("_t1-ce") != std::string::npos) || (fileName_wrap.find("_t1-gd") != std::string::npos) ||
    (fileName_wrap.find("_t1gd") != std::string::npos) || (fileName_wrap.find("t1gd_") != std::string::npos) ||
    (fileName_wrap.find("t1ce") != std::string::npos) || (fileName_wrap.find("t1-gad") != std::string::npos) ||
    (fileName_wrap.find("t1-ce") != std::string::npos) || (fileName_wrap.find("t1-gd") != std::string::npos) ||
    (fileName_wrap.find("t1ce.nii.gz") != std::string::npos) || (fileName_wrap.find("t1gd.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_T1CE;
  }
  else if ((fileName_wrap.find("_t1") != std::string::npos) || (fileName_wrap.find("t1_") != std::string::npos) || (fileName_wrap.find("t1.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_T1;
  }
  else if ((fileName_wrap.find("_t2") != std::string::npos) || (fileName_wrap.find("t2_") != std::string::npos) || (fileName_wrap.find("t2.nii.gz") != std::string::npos))
  {
    if ((fileName_wrap.find("flair") != std::string::npos))
    {
      ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR;
    }
    else
    {
      ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_T2;
    }
  }
  else if ((fileName_wrap.find("_flair") != std::string::npos) || (fileName_wrap.find("flair_") != std::string::npos) || (fileName_wrap.find(".flair.") != std::string::npos) 
    || (fileName_wrap.find("flair.nii.gz") != std::string::npos) || (fileName_wrap.find("fl.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_T2FLAIR;
  }
  else if ((fileName_wrap.find("_dti") != std::string::npos) || (fileName_wrap.find("dti_") != std::string::npos) ||
    (fileName_wrap.find("_b0") != std::string::npos) || (fileName_wrap.find("b0_") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_DTI;
  }
  else if ((fileName_wrap.find("_rad") != std::string::npos) || (fileName_wrap.find("rad_") != std::string::npos) ||
    (fileName_wrap.find("radial_") != std::string::npos) || (fileName_wrap.find("_radial") != std::string::npos) ||
    (fileName_wrap.find("radial.nii.gz") != std::string::npos) || (fileName_wrap.find("rad.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_RAD;
  }
  else if ((fileName_wrap.find("_ax") != std::string::npos) || (fileName_wrap.find("ax_") != std::string::npos) || 
    (fileName_wrap.find("axial_") != std::string::npos) || (fileName_wrap.find("_axial") != std::string::npos) ||
    (fileName_wrap.find("axial.nii.gz") != std::string::npos) || (fileName_wrap.find("ax.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_AX;
  }
  else if ((fileName_wrap.find("_fa") != std::string::npos) || (fileName_wrap.find("fa_") != std::string::npos) ||
    (fileName_wrap.find("fractional_") != std::string::npos) || (fileName_wrap.find("_fractional") != std::string::npos) ||
    (fileName_wrap.find("fractional.nii.gz") != std::string::npos) || (fileName_wrap.find("fa.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_FA;
  }
  else if ((fileName_wrap.find("_tr") != std::string::npos) || (fileName_wrap.find("tr_") != std::string::npos) ||
    (fileName_wrap.find("trace_") != std::string::npos) || (fileName_wrap.find("_trace") != std::string::npos) ||
    (fileName_wrap.find("trace.nii.gz") != std::string::npos) || (fileName_wrap.find("tr.nii.gz") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_TR;
  }
  else if ((fileName_wrap.find("_rcbv") != std::string::npos) || (fileName_wrap.find("rcbv_") != std::string::npos) || 
    (fileName_wrap.find("_rcbv_") != std::string::npos) || (fileName_wrap.find("rcbv") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_RCBV;
  }
  else if ((fileName_wrap.find("_psr") != std::string::npos) || (fileName_wrap.find("psr_") != std::string::npos) ||
    (fileName_wrap.find("_psr_") != std::string::npos) || (fileName_wrap.find("psr") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_PSR;
  }
  else if ((fileName_wrap.find("_ph") != std::string::npos) || (fileName_wrap.find("ph_") != std::string::npos) ||
    (fileName_wrap.find("_ph_") != std::string::npos) || (fileName_wrap.find("ph") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_PH;
  }
  else if ((fileName_wrap.find("perf") != std::string::npos) || (fileName_wrap.find("PERF") != std::string::npos) || 
    (fileName_wrap.find("DSC") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_PERFUSION;
  }
  else if ((fileName_wrap.find("label-map") != std::string::npos) || (fileName_wrap.find("label") != std::string::npos) ||
    (fileName_wrap.find("segmentation") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_SEG;
  }
  else if ((fileName_wrap.find("ct2pet") != std::string::npos) || (fileName_wrap.find("_ct.") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_CT;
  }
  else if ((fileName_wrap.find("_pet") != std::string::npos) || (fileName_wrap.find("pet_") != std::string::npos) || (fileName_wrap.find("_pet_") != std::string::npos))
  {
    ImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_PET;
  }
  return ImageSubType;
}

//! Checks for common file types that CaPTk can read properly, specially for DTI and DWI images
inline bool isExtensionSupported(const std::string inputExtension)
{
  if ((inputExtension == IMG_EXT || inputExtension == NII_EXT || inputExtension == NII_GZ_EXT))
  {
    return true;
  }
  else
  {
    return false;
  }
}

/** \brief Find total available memory (based on StackOverflow 2513505, Travis Gockel answer) */
inline unsigned long long getTotalInstalledMemory()
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

/** \brief Check if the total size of the files is more than a 
 * percentage of the available memory
 * */
inline bool isSizeOfLoadedFilesTooBig(QStringList files, std::string loggerFile = "", 
                               float maxPercentage = 0.05)
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
inline void WriteCSVFilesWithHorizontalAndVerticalHeaders(VariableSizeMatrixType inputdata, std::vector<std::string> horizontal_ids, std::vector<std::string> vertical_ids, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);

  myfile << " ";
  for(unsigned int col_index=0;col_index<vertical_ids.size();col_index++)
      myfile << ","<<vertical_ids[col_index];
  myfile << "\n";
  for (unsigned int row_index = 0; row_index < inputdata.Rows(); row_index++)
  {
   myfile << horizontal_ids[row_index];
    for (unsigned int col_index = 0; col_index < inputdata.Cols(); col_index++)
       myfile << "," << std::to_string(inputdata[row_index][col_index]);
    myfile << "\n";
  }
  myfile.close();
}


inline void WriteCSVFilesWithHorizontalAndVerticalHeaders(VariableSizeMatrixType inputdata, std::vector<std::string> horizontal_ids, std::string vertical_ids[], std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  myfile << " ";
  for(unsigned int col_index=0;col_index<vertical_ids->size();col_index++)
      myfile << ","<<vertical_ids[col_index];
  myfile << "\n";
  for (unsigned int row_index = 0; row_index < inputdata.Rows(); row_index++)
  {
   myfile << horizontal_ids[row_index];
    for (unsigned int col_index = 0; col_index < inputdata.Cols(); col_index++)
       myfile << "," << std::to_string(inputdata[row_index][col_index]);
    myfile << "\n";
  }
  myfile.close();
}

inline void WriteCSVFiles(VariableSizeMatrixType inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.Rows(); index1++)
  {
    for (unsigned int index2 = 0; index2 < inputdata.Cols(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(inputdata[index1][index2]);
      else
        myfile << "," << std::to_string(inputdata[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
}
inline void WriteCSVFiles(VectorVectorDouble inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.size(); index1++)
  {
    for (unsigned int index2 = 0; index2 < inputdata[0].size(); index2++)
    {
      if (index2 == 0)
        myfile << std::to_string(inputdata[index1][index2]);
      else
        myfile << "," << std::to_string(inputdata[index1][index2]);
    }
    myfile << "\n";
  }
  myfile.close();
}

inline void WriteCSVFiles(VariableLengthVectorType inputdata, std::string filepath, bool vertical=false)
{
  std::ofstream myfile;
  myfile.open(filepath);
  if (vertical == false)
  {
    for (unsigned int index1 = 0; index1 < inputdata.Size(); index1++)
      myfile << std::to_string(inputdata[index1]) << ",";

    myfile << "\n";
  }
  else
  {
    for (unsigned int index1 = 0; index1 < inputdata.Size(); index1++)
    {
      myfile << std::to_string(inputdata[index1]) << ",";
      if (index1 < inputdata.Size() - 1)
        myfile << "\n";
    }
  }
  myfile.close();
}
inline void WriteCSVFiles(std::vector<int> inputdata, std::string filepath)
{
  std::ofstream myfile;
  myfile.open(filepath);
  for (unsigned int index1 = 0; index1 < inputdata.size(); index1++)
    myfile << std::to_string(inputdata[index1]) << ",";

  myfile << "\n";
  myfile.close();
}
inline void WriteCSVFiles(std::vector<double> inputdata, std::string filepath,bool vertical=false)
{
  std::ofstream myfile;
  myfile.open(filepath);
  if (vertical == false)
  {
    for (unsigned int index1 = 0; index1 < inputdata.size(); index1++)
      myfile << std::to_string(inputdata[index1]) << ",";

    myfile << "\n";
  }
  else
  {
    for (unsigned int index1 = 0; index1 < inputdata.size(); index1++)
    {
      myfile << std::to_string(inputdata[index1]) << ",";
      if (index1 < inputdata.size() - 1)
        myfile << "\n";
    }
  }
  myfile.close();
}