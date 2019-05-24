/**
\file  cbicaCmdParser.cpp

\brief Implementation file for CmdParser class.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#if (_WIN32)
#define NOMINMAX
#include <direct.h>
#include <iostream>
#include <windows.h>
#include <conio.h>
#include <lmcons.h>
#include <Shlobj.h>
#include <filesystem>
#define GetCurrentDir _getcwd
//static bool WindowsDetected = true;
static const char  cSeparator = '\\';
//  static const char* cSeparators = "\\/";
#else
#include <dirent.h>
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <cstring>
#include <cstdlib>
#include <sys/types.h>
#include <errno.h>
#include <ftw.h>
#define GetCurrentDir getcwd
//static bool WindowsDetected = false;
static const char  cSeparator = '/';
//  static const char* cSeparators = "/";
#endif

#if (__APPLE__)
#include <mach-o/dyld.h>
#endif

#include <assert.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <vector>
#include "cbicaCmdParser.h"
#include "yaml-cpp/yaml.h"

#ifndef PROJECT_VERSION
#define PROJECT_VERSION "0.0.1"
#endif

/*
\namespace cbica
\brief Namespace for differentiating functions written for internal use
*/
namespace cbica
{
  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline bool directoryExists(const std::string &dName)
  {
    struct stat info;
    std::string dName_Wrap = dName;

    if (dName_Wrap[dName_Wrap.length() - 1] == '/')
    {
      dName_Wrap.erase(dName_Wrap.end() - 1);
    }

    if (stat(dName_Wrap.c_str(), &info) != 0)
      return false;
    else if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on windows
      return true;
    else
      return false;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline bool fileExists(const std::string &fName)
  {
    std::ifstream file_exists(fName.c_str());
    if (file_exists.good())
      return true;
    else
      return false;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline std::string getEnvironmentVariableValue(const std::string &environmentVariable)
  {
    std::string returnString = "";
    char tempValue[FILENAME_MAX];
#if defined(_WIN32)
    char tmp[FILENAME_MAX];
    size_t size = FILENAME_MAX;
    getenv_s(&size, tmp, size, environmentVariable.c_str()); // does not work, for some reason - needs to be tested
    std::string temp = cbica::stringReplace(tmp, "\\", "/");
    sprintf_s(tempValue, static_cast<size_t>(FILENAME_MAX), "%s", temp.c_str());
    tmp[0] = '\0';
#else
    char *tmp;
    tmp = std::getenv(environmentVariable.c_str());
    sprintf(tempValue, "%s", tmp);
#endif

    returnString = std::string(tempValue);
    tempValue[0] = '\0';

    return returnString;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline bool createDir(const std::string &dir_name)
  {
    //! Pure c++ based directory creation
#if defined(_WIN32)
    DWORD ftyp = GetFileAttributesA(dir_name.c_str()); // check if directory exists or not
    if (ftyp == INVALID_FILE_ATTRIBUTES)
      _mkdir(dir_name.c_str());
    return true;
#else
    DIR *pDir;
    pDir = opendir(dir_name.c_str()); // check if directory exists or not
    if (pDir == NULL)
      mkdir(dir_name.c_str(), 0777);
    return true;
#endif
    return false;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline char* constCharToChar(const std::string &input)
  {
    char *s = new char[input.size() + 1];
#ifdef _WIN32
    strcpy_s(s, input.size() + 1, input.c_str());
#else
    std::strcpy(s, input.c_str());
#endif
    return s;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline std::string iterateOverStringAndSeparators(const std::string &inputString, size_t &count, int enum_separator = 10)
  {
    std::string returnString = "";
    if (enum_separator == 10) // to get description
    {
      returnString = inputString.substr(count + 1); // get all characters after the separator was detected
    }
    else // for everything other than description
    {
      char testChar = inputString[count], separatorChar = *cbica::constCharToChar(getSeparator(enum_separator));
      size_t position, // position to start getting the substring
        separatorChecker = 2; // the configuration file needs the difference between two types of strings to be a single space (apart from the separator string)
      if (testChar == separatorChar)
      {
        count++;
        position = count;
        //stringStream.clear();
        //stringStream << inputString[count];
        //std::string testStr = stringStream.str(), testSep = getSeparator(enum_separator);
        //while (stringStream.str() != getSeparator(enum_separator))
        //{
        //  stringStream << inputString[count];
        //  returnString += stringStream.str();
        //  stringStream.clear();
        //  count++;
        //}
        testChar = inputString[count];
        while (testChar != separatorChar)
        {
          count++;
          testChar = inputString[count];
        }

        returnString = inputString.substr(position, count - position);
      }
      else // a small check as a contingency plan
      {
        while (separatorChecker > 0)
        {
          separatorChecker--;
          count++;
          testChar = inputString[count];
        }
      }
    }

    return returnString;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  std::string _getFullPath()
  {
#if defined(_WIN32)
    //! Initialize pointers to file and user names
    char path[FILENAME_MAX];
    GetModuleFileNameA(NULL, path, FILENAME_MAX);
    //_splitpath_s(filename, NULL, NULL, NULL, NULL, filename, NULL, NULL, NULL);
#elif __APPLE__
    char path[PATH_MAX];
    uint32_t size = PATH_MAX - 1;
    if (path != NULL)
    {
      if (_NSGetExecutablePath(path, &size) != 0)
      {
        std::cerr << "[getFullPath()] Error during getting full path..";
      }
    }
#else
    //! Initialize pointers to file and user names
    char path[PATH_MAX];
    if (::readlink("/proc/self/exe", path, sizeof(path) - 1) == -1)
      //path = dirname(path);
      std::cerr << "[getFullPath()] Error during getting full path..";
#endif

    std::string return_string = std::string(path);
    path[0] = '\0';

    return return_string;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline bool splitFileName(const std::string &dataFile, std::string &path,
    std::string &baseName, std::string &extension)
  {
    std::string dataFile_wrap = dataFile;
    std::vector< std::string > compressionFormats;
    compressionFormats.push_back(".gz");
    compressionFormats.push_back(".bz");
    compressionFormats.push_back(".zip");
    compressionFormats.push_back(".bz2");

    // check for compression formats
    for (size_t i = 0; i < compressionFormats.size(); i++)
    {
      if (dataFile_wrap.find(compressionFormats[i]) != std::string::npos)
      {
        dataFile_wrap = cbica::stringReplace(dataFile_wrap, compressionFormats[i], "");
        std::string tempExt;
        cbica::splitFileName(dataFile_wrap, path, baseName, tempExt);
        extension = tempExt + compressionFormats[i];
      }
    }
    if (!path.empty() && !baseName.empty() && !extension.empty())
    {
      return true;
    }
    else
    {
      //! Initialize pointers to file and user names
#if (_MSC_VER >= 1700)
      char basename_var[FILENAME_MAX], ext[FILENAME_MAX], path_name[FILENAME_MAX], drive_letter[FILENAME_MAX];
      //_splitpath(dataFile_wrap.c_str(), NULL, path_name, basename_var, ext);
      _splitpath_s(dataFile.c_str(), drive_letter, FILENAME_MAX, path_name, FILENAME_MAX, basename_var, FILENAME_MAX, ext, FILENAME_MAX);
#else
      char *basename_var, *ext, *path_name;
      path_name = dirname(cbica::constCharToChar(dataFile_wrap.c_str()));
      basename_var = basename(cbica::constCharToChar(dataFile_wrap.c_str()));
      ext = strrchr(cbica::constCharToChar(dataFile_wrap.c_str()), '.');
#endif

      //path sanity check
      if (path_name == NULL)
      {
        std::cerr << "No filename path has been detected.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        path =
#ifdef _WIN32
          std::string(drive_letter) +
#endif
          std::string(path_name);
      }
      path = cbica::stringReplace(path, "\\", "/"); // normalize path for Windows

      //base name sanity check
      if (basename_var == NULL)
      {
        std::cerr << "No filename base has been detected.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        baseName = std::string(basename_var);
      }

      //extension sanity check
      if (ext == NULL)
      {
        extension = "";
      }
      else
      {
        extension = std::string(ext);
      }

#if (_MSC_VER >= 1700)
      path_name[0] = NULL;
      basename_var[0] = NULL;
      ext[0] = NULL;
      drive_letter[0] = NULL;
#endif
      if (path[path.length() - 1] != '/')
      {
        path += "/";
      }

      return true;
    }
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  std::string _getExecutablePath()
  {
    std::string path, base, ext;
    cbica::splitFileName(cbica::_getFullPath(), path, base, ext);
    return path;
  }

  static inline std::string getFilenameBase(const std::string &filename, bool checkFile)
  {
    if (checkFile)
    {
      if (!fileExists(filename))
      {
        std::cerr << "[getFilenameBase()] Supplied file name'" << filename << "'wasn't found.\n";
        exit(EXIT_FAILURE);
      }
    }
    std::string path, base, ext;
    splitFileName(filename, path, base, ext);

    return base;
  }

  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  static inline std::string getExecutableName()
  {
    std::string return_string;
#if defined(_WIN32)
    //! Initialize pointers to file and user names
    char filename[FILENAME_MAX];
    GetModuleFileNameA(NULL, filename, FILENAME_MAX);
    return_string = getFilenameBase(filename, false);
    //_splitpath_s(filename, NULL, NULL, NULL, NULL, filename, NULL, NULL, NULL);
#elif __APPLE__
    char path[PATH_MAX];
    uint32_t size = PATH_MAX - 1;
    if (path != NULL)
    {
      if (_NSGetExecutablePath(path, &size) == 0)
      {
        return_string = getFilenameBase(std::string(path), false);
      }
    }
#else
    return_string = getEnvironmentVariableValue("_");
    return_string = stringReplace(return_string, "./", "");
    char path[PATH_MAX];
    if (path != NULL)
    {
      auto ret = ::readlink("/proc/self/exe", path, sizeof(path) - 1);
      if (ret == -1)
      {
        //free(path);
        path[0] = '\0';
      }
      path[ret] = 0;
    }
    auto temp = std::string(path);
    return_string = getFilenameBase(temp, false);
    path[0] = '\0';
#endif

    return return_string;
  }

  static inline std::string makeTempDir()
  {
    std::string returnDir = "", tempCheck, homeEnv;
#if defined(_WIN32)
    homeEnv = "USERPROFILE";
#else
    homeEnv = "HOME";
#endif    

    tempCheck = cbica::getEnvironmentVariableValue(homeEnv);
    tempCheck += "/tmp";

    if (cbica::directoryExists(tempCheck))
    {
      for (size_t i = 1; i <= FILENAME_MAX; i++)
      {
        returnDir = tempCheck + std::to_string(i);
        if (!cbica::directoryExists(returnDir))
        {
          break;
        }
      }
    }
    else
    {
      returnDir = tempCheck;
    }

    returnDir += "/";
    if (!createDir(returnDir))
    {
      std::cerr << "Could not create the temporary directory '" << returnDir << "'\n";
      exit(EXIT_FAILURE);
    }

    return returnDir;
  }


  void CmdParser::initializeClass(int &input_argc, std::vector< std::string > &input_argv, const std::string &input_exeName)
  {
#ifdef PROJECT_VERSION
    m_version = std::string(PROJECT_VERSION);
#else
    m_version = 0.1.0;
#endif    
    std::string path, base, ext;
    cbica::splitFileName(cbica::_getFullPath(), path, base, ext);
    if (input_exeName.empty())
    {
      m_exeName = getExecutableName();
    }
    else
    {
      m_exeName = input_exeName;
    }
    m_exePath = path;

    m_argc = input_argc;
    m_argv = input_argv;

    m_maxLength = 0;
    checkMaxLen = false;
    helpRequested = false;
    firstRun = true;
    argc1ignore = false;
    m_exampleOfUsage = "";

    m_optionalParameters.push_back(Parameter("u", "usage", cbica::Parameter::NONE, "", "Prints basic usage message", "", "", "", ""));
    m_optionalParameters.push_back(Parameter("h", "help", cbica::Parameter::NONE, "", "Prints verbose usage information", "", "", "", ""));
    m_optionalParameters.push_back(Parameter("v", "version", cbica::Parameter::NONE, "", "Prints information about software version", "", "", "", ""));
    m_optionalParameters.push_back(Parameter("rt", "run-test", cbica::Parameter::NONE, "", "Runs the tests", "", "", "", ""));
    m_optionalParameters.push_back(Parameter("cwl", "cwl", cbica::Parameter::NONE, "", "Generates a .cwl file for the software", "", "", "", ""));
  }

  CmdParser::CmdParser(int argc, char **argv, const std::string &exe_name)
  {
    if (m_argv.empty())
    {
      for (int i = 0; i < argc; i++)
      {
        m_argv.push_back(std::string(argv[i]));
      }
    }
    initializeClass(argc, m_argv, exe_name);
  }

  CmdParser::~CmdParser()
  {

  }

  static inline std::string getCurrentYear()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    localtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%d", timeinfo.tm_year + 1900);
#else
    tm *time_struct = localtime(&timer);
    sprintf(buffer, "%d", time_struct->tm_year + 1900);
#endif
    return std::string(buffer);
  }

  inline void copyrightNotice()
  {
    std::cout <<
      "\n==========================================================================\n" <<
      "Contact: software@cbica.upenn.edu\n\n" <<
      "Copyright (c) " << cbica::getCurrentYear() << " University of Pennsylvania. All rights reserved.\n" <<
      "See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html" <<
      "\n==========================================================================\n";
  }

  inline void CmdParser::getMaxLength()
  {
    m_minVerboseLength = 1024;
    m_maxVerboseLength = 0;
    m_maxLaconicLength = 0;
    m_maxLength = 0; // maximum length of laconic + verbose

    // loop through optional and required parameters separately
    for (size_t i = 0; i<m_optionalParameters.size(); ++i)
    {
      m_maxLength = m_maxLength < m_optionalParameters[i].length ? m_optionalParameters[i].length : m_maxLength;
      m_minVerboseLength = m_minVerboseLength > m_optionalParameters[i].verbose.length() ? m_optionalParameters[i].verbose.length() : m_minVerboseLength;
      m_maxVerboseLength = m_maxVerboseLength < m_optionalParameters[i].verbose.length() ? m_optionalParameters[i].verbose.length() : m_minVerboseLength;
      m_maxLaconicLength = m_maxLaconicLength < m_optionalParameters[i].laconic.length() ? m_optionalParameters[i].laconic.length() : m_maxLaconicLength;
    }

    for (size_t i = 0; i < m_requiredParameters.size(); ++i)
    {
      m_maxLength = m_maxLength < m_requiredParameters[i].length ? m_requiredParameters[i].length : m_maxLength;
      m_minVerboseLength = m_minVerboseLength > m_requiredParameters[i].verbose.length() ? m_requiredParameters[i].verbose.length() : m_minVerboseLength;
      m_maxVerboseLength = m_maxVerboseLength < m_requiredParameters[i].verbose.length() ? m_requiredParameters[i].verbose.length() : m_minVerboseLength;
      m_maxLaconicLength = m_maxLaconicLength < m_requiredParameters[i].laconic.length() ? m_requiredParameters[i].laconic.length() : m_maxLaconicLength;
    }

    m_maxLength += 5;

    checkMaxLen = true; // trigger flag for future checks

    if (!helpRequested && (m_argc != 1))
    {
      for (size_t i = 0; i<m_requiredParameters.size(); ++i)
      {
        // check if current required parameter has been supplied in the command line (obtained from argv)
        int tempPos;
        if (!CmdParser::compareParameter(m_requiredParameters[i].laconic, tempPos) && !helpRequested)
        {
          std::cout << "The required parameter '" << m_requiredParameters[i].laconic << "' is missing from the command line arguments you provided. See '" <<
            m_exeName << " --help' for extended help.\n\n";

          std::string m_exeName_wrap;
#ifdef _WIN32
          m_exeName_wrap = m_exeName + ".exe";
#else
          m_exeName_wrap = m_exeName;
#endif

          std::cout << "An exemplary usage scenario: \n\n" << m_exeName_wrap << " " << m_exampleOfUsage << "\n\n";

          exit(EXIT_FAILURE);
        }
      }
    }
  }

  void CmdParser::addOptionalParameter(const std::string &laconic, const std::string &verbose, const int &expectedDataType, const std::string &dataRange,
    const std::string &description_line1,
    const std::string &description_line2,
    const std::string &description_line3,
    const std::string &description_line4,
    const std::string &description_line5)
  {
    if ((laconic == "u") || (laconic == "h") || (laconic == "v") || (laconic == "rt") || (laconic == "cwl"))
    {
      return;
    }
    if (laconic == "")
    {
      std::cerr << "Laconic parameter cannot be empty";
      exit(EXIT_FAILURE);
    }
    if (verbose == "")
    {
      std::cerr << "Verbose parameter cannot be empty";
      exit(EXIT_FAILURE);
    }
    if (description_line1 == "")
    {
      std::cerr << "Failure to initialize an empty string as description_line1";
      exit(EXIT_FAILURE);
    }

    std::string verbose_wrap = verbose;
    stringReplace(verbose_wrap, " ", "");

    m_optionalParameters.push_back(Parameter(laconic, verbose_wrap, expectedDataType, dataRange, description_line1, description_line2, description_line3, description_line4, description_line5));
  }

  void CmdParser::addRequiredParameter(const std::string &laconic, const std::string &verbose, const int &expectedDataType, const std::string &dataRange,
    const std::string &description_line1,
    const std::string &description_line2,
    const std::string &description_line3,
    const std::string &description_line4,
    const std::string &description_line5)
  {
    //std::cout << laconic << verbose << verbose << dataRange << description_line1 << std::endl;
    if ((laconic == "u") || (laconic == "h") || (laconic == "v") || (laconic == "rt") || (laconic == "cwl"))
    {
      return;
    }
    if (laconic == "")
    {
      std::cerr << "Laconic parameter cannot be empty";
      exit(EXIT_FAILURE);
    }
    if (verbose == "")
    {
      std::cerr << "Verbose parameter cannot be empty";
      exit(EXIT_FAILURE);
    }
    if (description_line1 == "")
    {
      std::cerr << "Failure to initialize an empty string as description_line1";
      exit(EXIT_FAILURE);
    }

    std::string verbose_wrap = verbose;
    stringReplace(verbose_wrap, " ", "");

    m_requiredParameters.push_back(Parameter(laconic, verbose_wrap, expectedDataType, dataRange, description_line1, description_line2, description_line3, description_line4, description_line5));
  }

  void CmdParser::addParameter(const std::string &laconic, const std::string &verbose, const int &expectedDataType, const std::string &dataRange,
    const std::string &description_line1,
    const std::string &description_line2,
    const std::string &description_line3,
    const std::string &description_line4,
    const std::string &description_line5)
  {
    CmdParser::addOptionalParameter(laconic, verbose, expectedDataType, dataRange, description_line1, description_line2, description_line3, description_line4, description_line5);
  }

  inline void CmdParser::writeParameters(const std::vector< Parameter > &inputParameters, bool verbose)
  {
    std::string spaces_verb_line2;
    for (size_t n = 0; n < m_maxLength + 6; n++)
    {
      spaces_verb_line2.append(" ");
    }


    for (size_t i = 0; i < inputParameters.size(); i++)
    {
      std::string spaces_lac, spaces_verb;

      for (size_t n = 0; n < m_maxLaconicLength - inputParameters[i].laconic.length(); n++)
      {
        spaces_lac.append(" ");
      }

      for (size_t n = 0; n < m_maxLength - inputParameters[i].verbose.length() - spaces_lac.length() - 3; n++)
      {
        spaces_verb.append(" ");
      }
      if (inputParameters[i].laconic.length() == 1)
      {
        spaces_verb.append(" ");
      }

      if (spaces_lac.empty() && !spaces_verb.empty())
      {
        for (size_t n = 0; n <= m_maxLaconicLength - inputParameters[i].laconic.length(); n++)
        {
          spaces_verb.pop_back();
        }
      }

      std::cout << "[" << spaces_lac << "-" << inputParameters[i].laconic << ", --" <<
        inputParameters[i].verbose << "]" << spaces_verb <<
        inputParameters[i].descriptionLine1 << " \n";

      if (inputParameters[i].descriptionLine2 != "")
      {
        std::cout << spaces_verb_line2 << inputParameters[i].descriptionLine2 << "\n";
        if (inputParameters[i].descriptionLine3 != "")
        {
          std::cout << spaces_verb_line2 << inputParameters[i].descriptionLine3 << "\n";
          if (inputParameters[i].descriptionLine4 != "")
          {
            std::cout << spaces_verb_line2 << inputParameters[i].descriptionLine4 << "\n";
            if (inputParameters[i].descriptionLine5 != "")
            {
              std::cout << spaces_verb_line2 << inputParameters[i].descriptionLine5 << "\n";
            }
          }
        }
      }

      if (verbose && (inputParameters[i].laconic != "u") && (inputParameters[i].laconic != "h") && (inputParameters[i].laconic != "v") && (inputParameters[i].laconic != "rt") && (inputParameters[i].laconic != "cwl") )
      {
        std::cout << spaces_verb_line2 << "Expected Type  :: " << inputParameters[i].dataType_string << "\n" <<
          spaces_verb_line2 << "Expected Range :: " << inputParameters[i].dataRange << "\n";
      }

      std::cout << "\n"; // an extra to keep coherence
    }
  }

  void CmdParser::echoUsage()
  {
    if (!checkMaxLen)
    {
      getMaxLength();
    }
    std::cout << "Executable Name: " << m_exeName << " v" << m_version
      << "\n\n" << "Usage:\n\n";

    std::cout << "Required parameters:\n\n";

    writeParameters(m_requiredParameters, false);
    std::cout << "Optional parameters:\n\n";
    writeParameters(m_optionalParameters, false);

    copyrightNotice();
  }

  void CmdParser::echoHelp()
  {
    if (!checkMaxLen)
    {
      getMaxLength();
    }
    std::cout << "Executable Name: " << m_exeName << " v" << m_version
      << "\n\n";
    
    std::cout << "Description:\n\n" << m_description << "\n\n";

    std::cout << "Usage:\n\n";
    std::cout << ":::Required parameters:::\n\n";
    writeParameters(m_requiredParameters, true);
    std::cout << ":::Optional parameters:::\n\n";
    writeParameters(m_optionalParameters, true);

    std::string m_exeName_wrap;
#ifdef _WIN32
    m_exeName_wrap = m_exeName + ".exe";
#else
    m_exeName_wrap = m_exeName;
#endif

    if (m_exampleOfUsage != "") // this is deprecated
    {
      std::cout << "For example: \n\n" <<
        m_exeName_wrap << " " << m_exampleOfUsage << "\n";
    }

    if (!m_exampleUsageAndDescription.empty()) // this will become required
    {
      std::cout << "Examples of Usage:\n\n";

      for (size_t i = 0; i < m_exampleUsageAndDescription.size(); i++)
      {
        std::cout << " Command: " << m_exeName_wrap << " " << m_exampleUsageAndDescription[i].first << "\n";
        std::cout << " Result : " << m_exampleUsageAndDescription[i].second << "\n\n";
      }
      std::cout << "\n";
    }

    copyrightNotice();
  }

  void CmdParser::echoVersion()
  {
    std::cout << "Executable Name: " << m_exeName << "\n" << "        Version: " <<
      m_version << "\n";

    copyrightNotice();
  }

  inline std::string internal_compare(const std::string &check_string, const int check_length)
  {
    switch (std::abs(static_cast<int>(check_string.length() - check_length)))
    {
    case 1:
      return ("-" + check_string);
      break;
    case 2:
      return ("--" + check_string);
      break;
    default:
      return (check_string);
      break;
    }
  }

  inline void CmdParser::verbose_check(std::string &input_string)
  {
    std::string input_string_lower = input_string;
    std::transform(input_string_lower.begin(), input_string_lower.end(), input_string_lower.begin(), ::tolower);
    if ((input_string_lower == "usage") || (input_string_lower == "-usage") || (input_string_lower == "--usage")
      || (input_string_lower == "u") || (input_string_lower == "-u") || (input_string_lower == "--u"))
    {
      input_string = "u";
    }
    else if ((input_string_lower == "help") || (input_string_lower == "-help") || (input_string_lower == "--help")
      || (input_string_lower == "h") || (input_string_lower == "-h") || (input_string_lower == "--h"))
    {
      input_string = "h";
    }
    else if ((input_string_lower == "version") || (input_string_lower == "-version") || (input_string_lower == "--version")
      || (input_string_lower == "v") || (input_string_lower == "-v") || (input_string_lower == "--v"))
    {
      input_string = "v";
    }
    else if ((input_string_lower == "run-test") || (input_string_lower == "-run-test") || (input_string_lower == "--run-test")
      || (input_string_lower == "rt") || (input_string_lower == "-rt") || (input_string_lower == "--rt"))
    {
      input_string = "rt";
    }
    else if ((input_string_lower == "cwl") || (input_string_lower == "-cwl") || (input_string_lower == "--cwl"))
    {
      input_string = "cwl";
    }


    if (!checkMaxLen)
    {
      getMaxLength();
    }
    if (input_string.length() > m_maxLaconicLength)
    {
      input_string = cbica::stringReplace(input_string, "--", "");
      input_string = cbica::stringReplace(input_string, "-", "");

      for (size_t i = 0; i < m_requiredParameters.size(); i++)
      {
        input_string = m_requiredParameters[i].verbose == input_string ? m_requiredParameters[i].laconic : input_string;
      }

      for (size_t i = 0; i < m_optionalParameters.size(); i++)
      {
        input_string = m_optionalParameters[i].verbose == input_string ? m_optionalParameters[i].laconic : input_string;
      }

      return;
    }
  }

  bool CmdParser::compareParameter(const std::string &execParamToCheck, int &position, bool automaticEcho)
  {
    // check for argc values during the first run otherwise don't
    if (firstRun)
    {
      //assert(("Example of usage and respective description is needed; use 'parser.addExampleUsage()'", !m_exampleUsageAndDescription.empty());
      //assert(("Application description is needed; use 'parser.addApplicationDescription()'", !m_description.empty());
      if (m_argc > static_cast< int >(2 * (m_optionalParameters.size() + m_requiredParameters.size() - 3) + 1))
      {
        std::cerr << "Extra parameters passed, please check usage. Exiting.\n\n";
        echoUsage();
        exit(EXIT_FAILURE);
      }

      if (!argc1ignore)
      {
        if (m_argc < 2)
        {
          std::cerr << "Insufficient parameters passed, please check usage. Exiting.\n\n";
          echoUsage();
          exit(EXIT_FAILURE);
        }
      }

      firstRun = false;
    }

    std::string execParamToCheck_wrap = execParamToCheck;
    verbose_check(execParamToCheck_wrap);

    for (int i = 1; i < m_argc; i++)
    {
      std::string inputParamToCheck = m_argv[i];
      verbose_check(inputParamToCheck);
      if (inputParamToCheck == "u")
      {
        helpRequested = true;
        position = i;
        if (automaticEcho)
        {
          echoUsage();
          exit(EXIT_SUCCESS);
        }
        return true;
      }
      if (inputParamToCheck == "h")
      {
        helpRequested = true;
        position = i;
        if (automaticEcho)
        {
          echoHelp();
          exit(EXIT_SUCCESS);
        }
        return true;
      }
      if (inputParamToCheck == "v")
      {
        helpRequested = true;
        position = i;
        if (automaticEcho)
        {
          echoVersion();
          exit(EXIT_SUCCESS);
        }
        return true;
      }
      if (inputParamToCheck == "rt")
      {
        helpRequested = true;
        position = i;
        // writeCWLFile(_getExecutablePath(), false);
        exit(EXIT_SUCCESS);
        //return true;
      }
      if (inputParamToCheck == "cwl")
      {
        helpRequested = true;
        position = i;
        //std::cout << "[DEBUG] m_exePath: "<< m_exePath << "\n";
        writeCWLFile(m_exePath, true);
        exit(EXIT_SUCCESS);
        return true;
      }
      if (!checkMaxLen)
      {
        getMaxLength();
      }

      if (inputParamToCheck == execParamToCheck_wrap)
      {
        position = i;
        return true;
      }
      else
      {
        std::string inputCheck, execCheck;
        const unsigned int minLength = static_cast<unsigned int>(std::max(
          inputParamToCheck.length(), execParamToCheck_wrap.length()));

        inputCheck = internal_compare(inputParamToCheck, minLength);
        execCheck = internal_compare(execParamToCheck_wrap, minLength);

        if (inputCheck == execCheck)
        {
          position = i;
          return true;
        }
      }
    }

    return false;
  }

  bool CmdParser::compareParameter(const std::string &execParamToCheck, bool automaticEcho)
  {
    int position;
    return compareParameter(execParamToCheck, position, automaticEcho);
  }

  bool CmdParser::isPresent(const std::string &execParamToCheck, bool automaticEcho)
  {
    return compareParameter(execParamToCheck, automaticEcho);
  }

  std::string CmdParser::getDescription(const std::string &execParamToCheck, bool NewLine = false)
  {
    int noMoreChecks = 0; // ensures that extra checks are not done for parameters
    if (execParamToCheck == "")
    {
      std::cerr << "Parameter cannot be an empty string. Please try again.\n";
      exit(EXIT_FAILURE);
    }
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    size_t i = 0;
    while ((i < m_requiredParameters.size()) && (noMoreChecks < 1))
    {
      if ((m_requiredParameters[i].laconic == execParamToCheck) ||
        (m_requiredParameters[i].verbose == execParamToCheck))
      {
        if (NewLine)
        {
          return m_requiredParameters[i].descriptionLine1 + "\n" + m_requiredParameters[i].descriptionLine2 + "\n" +
            m_requiredParameters[i].descriptionLine3 + "\n" + m_requiredParameters[i].descriptionLine4 + "\n" +
            m_requiredParameters[i].descriptionLine5;
        }
        else
        {
          return m_requiredParameters[i].descriptionLine1 + " " + m_requiredParameters[i].descriptionLine2 + " " +
            m_requiredParameters[i].descriptionLine3 + " " + m_requiredParameters[i].descriptionLine4 + " " +
            m_requiredParameters[i].descriptionLine5;
        }
      }
      i++;
    }

    i = 0;
    while ((i < m_optionalParameters.size()) && (noMoreChecks < 1))
    {
      if ((m_optionalParameters[i].laconic == execParamToCheck) ||
        (m_optionalParameters[i].verbose == execParamToCheck))
      {
        if (NewLine)
        {
          return m_optionalParameters[i].descriptionLine1 + "\n" + m_optionalParameters[i].descriptionLine2 + "\n" +
            m_optionalParameters[i].descriptionLine3 + "\n" + m_optionalParameters[i].descriptionLine4 + "\n" +
            m_optionalParameters[i].descriptionLine5;
        }
        else
        {
          return m_optionalParameters[i].descriptionLine1 + " " + m_optionalParameters[i].descriptionLine2 + " " +
            m_optionalParameters[i].descriptionLine3 + " " + m_optionalParameters[i].descriptionLine4 + " " +
            m_optionalParameters[i].descriptionLine5;
        }
        noMoreChecks = 1;
      }
      i++;
    }

    return "";
  }

  std::string CmdParser::getDataTypeAsString(const std::string &execParamToCheck)
  {
    int noMoreChecks = 0; // ensures that extra checks are not done for parameters
    if (execParamToCheck == "")
    {
      std::cerr << "Parameter cannot be an empty string. Please try again.\n";
      exit(EXIT_FAILURE);
    }
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    size_t i = 0;
    while ((i < m_requiredParameters.size()) && (noMoreChecks < 1))
    {
      if ((m_requiredParameters[i].laconic == execParamToCheck) ||
        (m_requiredParameters[i].verbose == execParamToCheck))
      {
        return m_requiredParameters[i].dataType_string;
        noMoreChecks = 1;
      }
      i++;
    }

    i = 0;
    while ((i < m_optionalParameters.size()) && (noMoreChecks < 1))
    {
      if ((m_optionalParameters[i].laconic == execParamToCheck) ||
        (m_optionalParameters[i].verbose == execParamToCheck))
      {
        return m_optionalParameters[i].dataType_string;
        noMoreChecks = 1;
      }
      i++;
    }

    return "";
  }

  int CmdParser::getDataTypeAsEnumCode(const std::string &execParamToCheck)
  {
    bool noMoreChecks = false; // ensures that extra checks are not done for parameters
    if (execParamToCheck == "")
    {
      std::cerr << "Parameter cannot be an empty string. Please try again.\n";
      exit(EXIT_FAILURE);
    }
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    size_t i = 0;
    while ((i < m_requiredParameters.size()) && !noMoreChecks)
    {
      if ((m_requiredParameters[i].laconic == execParamToCheck) ||
        (m_requiredParameters[i].verbose == execParamToCheck))
      {
        return m_requiredParameters[i].dataType_enumCode;
        noMoreChecks = true;
      }
      i++;
    }

    i = 0;
    while ((i < m_optionalParameters.size()) && !noMoreChecks)
    {
      if ((m_optionalParameters[i].laconic == execParamToCheck) ||
        (m_optionalParameters[i].verbose == execParamToCheck))
      {
        return m_optionalParameters[i].dataType_enumCode;
        noMoreChecks = true;
      }
      i++;
    }

    return -1;
  }

  void CmdParser::getParameterValue(const std::string &execParamToCheck, bool &parameterValue)
  {
    if (getDataTypeAsEnumCode(execParamToCheck) != cbica::Parameter::BOOLEAN)
    {
      std::cerr << "The data type of the requested parameter, '" << execParamToCheck << "' is classified as '" << getDataTypeAsString(execParamToCheck) <<
        "' and cannot be returned as a BOOL.\n";
      exit(EXIT_FAILURE);
    }
    int position;
    if (compareParameter(execParamToCheck, position))
    {
      if (position < (m_argc - 1))
      {
        std::string rawValue = m_argv[position + 1];
        if ((rawValue == "1") || (rawValue == "true") || (rawValue == "True") || (rawValue == "TRUE") ||
          (rawValue == "yes") || (rawValue == "Yes") || (rawValue == "YES") ||
          (rawValue.empty()) || (rawValue[0] == '-')) // if the parameter is just passed as a flag, assume that the user wants it enabled; '-' check is basically to check if the next parameter starts
        {
          parameterValue = true; // return value is a bool
          return;
        }
        else
        {
          parameterValue = false; // return value is a bool
          return;
        }
      }
      else
      {
        parameterValue = true; // the parameter has been passed as a flag at the end of the command
        return;
      }
    }
  }

  void CmdParser::getParameterValue(const std::string &execParamToCheck, int &parameterValue)
  {
    if (getDataTypeAsEnumCode(execParamToCheck) != cbica::Parameter::INTEGER)
    {
      std::cerr << "The data type of the requested parameter, '" << execParamToCheck << "' is classified as '" << getDataTypeAsString(execParamToCheck) <<
        "' and cannot be returned as a INTEGER.\n";
      exit(EXIT_FAILURE);
    }
    int position;
    if (compareParameter(execParamToCheck, position))
    {
      parameterValue = std::atoi(m_argv[position + 1].c_str()); // return value is an integer
      return;
    }
    else
    {
      parameterValue = -1;
      return;
    }
  }

  void CmdParser::getParameterValue(const std::string &execParamToCheck, size_t &parameterValue)
  {
    if (getDataTypeAsEnumCode(execParamToCheck) != cbica::Parameter::INTEGER)
    {
      std::cerr << "The data type of the requested parameter, '" << execParamToCheck << "' is classified as '" << getDataTypeAsString(execParamToCheck) <<
        "' and cannot be returned as a INTEGER.\n";
      exit(EXIT_FAILURE);
    }
    int position;
    if (compareParameter(execParamToCheck, position))
    {
      parameterValue = std::atoi(m_argv[position + 1].c_str()); // return value is an integer
      return;
    }
    else
    {
      parameterValue = 0;
      return;
    }
  }

  void CmdParser::getParameterValue(const std::string &execParamToCheck, float &parameterValue)
  {
    if (getDataTypeAsEnumCode(execParamToCheck) != cbica::Parameter::FLOAT)
    {
      std::cerr << "The data type of the requested parameter, '" << execParamToCheck << "' is classified as '" << getDataTypeAsString(execParamToCheck) <<
        "' and cannot be returned as a FLOAT.\n";
      exit(EXIT_FAILURE);
    }
    int position;
    if (compareParameter(execParamToCheck, position))
    {
      parameterValue = static_cast<float>(std::atof(m_argv[position + 1].c_str())); // return value is a float
      return;
    }
    else
    {
      parameterValue = -1;
      return;
    }
  }

  void CmdParser::getParameterValue(const std::string &execParamToCheck, std::string &parameterValue)
  {
    int returnCode = getDataTypeAsEnumCode(execParamToCheck);
    if ((returnCode != cbica::Parameter::STRING))
    {
      if (!((returnCode == cbica::Parameter::NONE) || // check if type is NONE or FILE or DIR, if yes then it is not an error
        (returnCode == cbica::Parameter::FILE) ||
        (returnCode == cbica::Parameter::DIRECTORY)))
      {
        std::cerr << "The data type of the requested parameter, '" << execParamToCheck << "' is classified as '" << getDataTypeAsString(execParamToCheck) <<
          "' and cannot be returned as a STRING.\n";
        exit(EXIT_FAILURE);
      }
    }
    int position;
    if (compareParameter(execParamToCheck, position))
    {
      parameterValue = m_argv[position + 1]; // return value is a string
      return;
    }
    else
    {
      parameterValue = "";
      return;
    }
  }

  void CmdParser::exampleUsage(const std::string &usageOfExe)
  {
    std::cerr << "!!! CmdParser::exampleUsage() IS DEPRECATED; use CmdParser::AddExampleUsage instead !!!\n";
    m_exampleOfUsage = usageOfExe;
    m_exampleOfUsage = cbica::stringReplace(m_exampleOfUsage, m_exeName + ".exe", "");
    m_exampleOfUsage = cbica::stringReplace(m_exampleOfUsage, "./" + m_exeName, "");
  }

  void CmdParser::addExampleUsage(const std::string &commandExcludingExeName, const std::string &descriptionOfCommand)
  {
    auto tempUsage = commandExcludingExeName;

    // remove Windows/Linux/macOS specific usage call pattern
    tempUsage = cbica::stringReplace(tempUsage, m_exeName + ".exe", "");
    tempUsage = cbica::stringReplace(tempUsage, "./" + m_exeName, "");
    tempUsage = cbica::stringReplace(tempUsage, "/" + m_exeName, "");

    m_exampleUsageAndDescription.push_back(std::make_pair(tempUsage, descriptionOfCommand));
  }

  void CmdParser::addApplicationDescription(const std::string &description)
  {
    m_description = description;
  }

  void CmdParser::writeCWLFile(const std::string &dirName, bool overwriteFile = true) 
  {
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    std::string dirName_wrap;
    if (!cbica::directoryExists(dirName) || (dirName.empty()))
    {
      dirName_wrap = cbica::makeTempDir();
    }
    dirName_wrap = cbica::stringReplace(dirName, "\\", "/");
    if (dirName_wrap.substr(dirName_wrap.length() - 1) != "/")
    {
      dirName_wrap += "/";
    }

    std::string cwlfileName = dirName_wrap + m_exeName + ".cwl";

    // std::cout << "[DEBUG]dirName_wrap: " << dirName_wrap << std::endl;
    // std::cout << "[DEBUG]m_exeName: " << m_exeName << std::endl;
    // std::cout << "[DEBUG]cwlfileName: " << cwlfileName << std::endl;
    
    std::ofstream file;
    if (!cbica::fileExists(cwlfileName) || overwriteFile)
    {
      file.open(cwlfileName.c_str());

      YAML::Node config = YAML::LoadFile(cwlfileName);

      config["cwlVersion"] = "v1.0";
      config["class"] = "CommandLineTool";
      config["version"] = m_version;
      config["baseCommand"] = (m_exeName);
      
      YAML::Node inputs = config["inputs"];
      
      for (size_t i = 0; i < m_requiredParameters.size(); i++)
      {
        config["inputs"]["-" + m_requiredParameters[i].verbose];
        config["inputs"][m_requiredParameters[i].verbose]["type"] =
          (m_requiredParameters[i].dataType_string == "STRING") ? "string" :
          (m_requiredParameters[i].dataType_string == "DIRECTORY") ? "Directory" :
          (m_requiredParameters[i].dataType_string == "FLOAT") ? "float" :
          (m_requiredParameters[i].dataType_string == "BOOL") ? "boolean" :
          (m_requiredParameters[i].dataType_string == "NONE") ? "string" :
          (m_requiredParameters[i].dataType_string == "UNKNOWN") ? "string" :
          (m_requiredParameters[i].dataType_string == "INTEGER") ? "int" :
          (m_requiredParameters[i].dataType_string == "FILE") ? "File" :
          "string";
        config["inputs"][m_requiredParameters[i].verbose]["label"] = m_requiredParameters[i].dataRange == "" ? "none" : m_requiredParameters[i].dataRange;
        YAML::Node inputBinding = config["inputBinding"];
        config["inputs"][m_requiredParameters[i].verbose]["inputBinding"]["position"] = 1;
        config["inputs"][m_requiredParameters[i].verbose]["inputBinding"]["prefix"] = "-" + m_requiredParameters[i].laconic;
        config["inputs"][m_requiredParameters[i].verbose]["doc"] =
          (m_requiredParameters[i].descriptionLine1 == "" ? "" : (m_requiredParameters[i].descriptionLine1 + ".")) +
          (m_requiredParameters[i].descriptionLine2 == "" ? "" : (m_requiredParameters[i].descriptionLine2 + ".")) +
          (m_requiredParameters[i].descriptionLine3 == "" ? "" : (m_requiredParameters[i].descriptionLine3 + ".")) +
          (m_requiredParameters[i].descriptionLine4 == "" ? "" : (m_requiredParameters[i].descriptionLine4 + ".")) +
          (m_requiredParameters[i].descriptionLine5 == "" ? "" : (m_requiredParameters[i].descriptionLine5 + "."));
      }

      if (m_optionalParameters.size() > 0) {
        for (size_t i = 0; i < m_optionalParameters.size(); i++)
        {
          if (m_optionalParameters[i].verbose == "help" ||
            m_optionalParameters[i].verbose == "usage" ||
            m_optionalParameters[i].verbose == "version" ||
            m_optionalParameters[i].verbose == "LogFile") {
            continue;
          }
          else {
            config["inputs"]["-" + m_optionalParameters[i].verbose];
            config["inputs"][m_optionalParameters[i].verbose]["type"] =
              (m_optionalParameters[i].dataType_string == "STRING") ? "string?" :
              (m_optionalParameters[i].dataType_string == "DIRECTORY") ? "Directory?" :
              (m_optionalParameters[i].dataType_string == "FLOAT") ? "float?" :
              (m_optionalParameters[i].dataType_string == "BOOL") ? "boolean?" :
              (m_optionalParameters[i].dataType_string == "NONE") ? "string?" :
              (m_optionalParameters[i].dataType_string == "UNKNOWN") ? "string?" :
              (m_optionalParameters[i].dataType_string == "INTEGER") ? "int?" :
              (m_optionalParameters[i].dataType_string == "FILE") ? "File?" :
              "string?";
            config["inputs"][m_optionalParameters[i].verbose]["label"] = m_optionalParameters[i].dataRange == "" ? "none" : m_optionalParameters[i].dataRange;
            config["inputs"][m_optionalParameters[i].verbose]["inputBinding"];
            config["inputs"][m_optionalParameters[i].verbose]["inputBinding"]["position"] = 1;
            config["inputs"][m_optionalParameters[i].verbose]["inputBinding"]["prefix"] = "-" + m_optionalParameters[i].laconic;
            config["inputs"][m_optionalParameters[i].verbose]["doc"] =
              (m_optionalParameters[i].descriptionLine1 == "" ? "" : (m_optionalParameters[i].descriptionLine1 + ".")) +
              (m_optionalParameters[i].descriptionLine2 == "" ? "" : (m_optionalParameters[i].descriptionLine2 + ".")) +
              (m_optionalParameters[i].descriptionLine3 == "" ? "" : (m_optionalParameters[i].descriptionLine3 + ".")) +
              (m_optionalParameters[i].descriptionLine4 == "" ? "" : (m_optionalParameters[i].descriptionLine4 + ".")) +
              (m_optionalParameters[i].descriptionLine5 == "" ? "" : (m_optionalParameters[i].descriptionLine5 + "."));
          }
        }
      }

      std::ofstream fout(cwlfileName);
      fout << config;
    }

    return;
  }

  std::string CmdParser::GetCommandFromCWL(const std::string &inpDir, const std::string & cwlDir)
  {
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    std::string inpDirName_warp, cwlDirName_warp;


    if (!cbica::directoryExists(inpDir) || (inpDir.empty()))
    {
      inpDirName_warp = cbica::makeTempDir();
    }
    inpDirName_warp = cbica::stringReplace(inpDir, "\\", "/");
    if (inpDirName_warp.substr(inpDirName_warp.length() - 1) != "/")
    {
      inpDirName_warp += "/";
    }

    //cwl dir
    if (!cbica::directoryExists(cwlDir) || (cwlDir.empty()))
    {
      cwlDirName_warp = cbica::makeTempDir();
    }
    cwlDirName_warp = cbica::stringReplace(cwlDir, "\\", "/");
    if (cwlDirName_warp.substr(cwlDirName_warp.length() - 1) != "/")
    {
      cwlDirName_warp += "/";
    }

    std::string inpFileName = inpDirName_warp + m_exeName + ".yml";
    std::string cwlFileName = cwlDirName_warp + m_exeName + ".cwl";

    YAML::Node input = YAML::LoadFile(inpFileName);
    YAML::Node config = YAML::LoadFile(cwlFileName);

    std::string cmd = m_exeName;

    for (YAML::const_iterator it = input.begin(); it != input.end(); ++it) {
      std::string key = it->first.as<std::string>();
      std::string value = it->second.as<std::string>();

      std::string laconic = getLaconic(key);
      if (laconic == "")
      {
        std::cerr << "No matching parameter found for '" << key << "'.\n";
        exit(EXIT_FAILURE);
      }

      if (!input[key]) {
        if (config["inputs"][key]["default"]) {
          value = config["inputs"][key]["default"].as<std::string>();
        }
        else {
          std::cerr << "There is no value specified for parameter'" << key << "'Default is also not specified.\n";
          exit(EXIT_FAILURE);
        }
      }
      cmd += " -" + laconic + " " + value;
    }

    logCWL(inpFileName, cwlFileName);
    return cmd;
  }

  void CmdParser::readCWLFile(const std::string &path_to_CWL_File, bool getDescription)
  {
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    std::string dirName_wrap;
    if (!cbica::directoryExists(path_to_CWL_File) || (path_to_CWL_File.empty()))
    {
      dirName_wrap = cbica::makeTempDir();
    }
    dirName_wrap = cbica::stringReplace(path_to_CWL_File, "\\", "/");
    if (dirName_wrap.substr(dirName_wrap.length() - 1) != "/")
    {
      dirName_wrap += "/";
    }

    std::string cwlfileName = dirName_wrap + m_exeName + ".cwl";

    YAML::Node config = YAML::LoadFile(cwlfileName);

    std::vector<Parameter> returnVector;

    if (!config)
    {
      std::cerr << "File '" << path_to_CWL_File << "' not found.\n";
      exit(EXIT_FAILURE);
    }
    std::string line;

    std::string parameter, parameterDataType, parameterDataRange, parameterDescription = "";

    YAML::Node characterType = config["inputs"];

    bool optional = false;
    int con = 0;
    for (YAML::const_iterator it = characterType.begin(); it != characterType.end(); ++it) {
      std::string key = it->first.as<std::string>();
      std::string verbose = key;

      std::string laconic = config["inputs"][key]["inputBinding"]["prefix"].as<std::string>();

      laconic.erase(std::remove(laconic.begin(), laconic.end(), '-'), laconic.end());


      std::string dataRange = config["inputs"][key]["label"].as<std::string>();
      std::string desc1 = config["inputs"][key]["doc"].as<std::string>();

      std::string dataType = config["inputs"][key]["type"].as<std::string>();

      if (dataType.find('?') != std::string::npos)
        optional = true;// find
      else
        optional = false;// not find

      cbica::Parameter::Type type;

      if (dataType == "File") {
        type = cbica::Parameter::FILE;
      }
      else if (dataType == "Directory")
      {
        type = cbica::Parameter::DIRECTORY;
      }
      else if (dataType == "string")
      {
        type = cbica::Parameter::STRING;
      }
      else if (dataType == "int")
      {
        type = cbica::Parameter::INTEGER;
      }
      else if (dataType == "float")
      {
        type = cbica::Parameter::FLOAT;
      }
      else if (dataType == "boolean")
      {
        type = cbica::Parameter::BOOLEAN;
      }
      else 
      {
        type = cbica::Parameter::NONE;
      }
      
      if (optional)
      {
        addOptionalParameter(laconic, verbose, type, dataRange, desc1);
      }
      else 
      {
        addRequiredParameter(laconic, verbose, type, dataRange, desc1);
      }

      con++;
    }
  }

  void CmdParser::createNode(const std::string & nodeString)
  {
    YAML::Node rootNode;
    if (!rootNode)
    {
      std::cerr << "Root node is empty\n";
      exit(EXIT_FAILURE);
    }

    YAML::Node node = rootNode[nodeString];

    return;
  }

  void CmdParser::deleteNode(const std::string & nodeString)
  {
    YAML::Node parent;
    if (!parent)
    {
      std::cerr << "Root node is empty\n";
      exit(EXIT_FAILURE);
    }

    //parent[nodeString].remove;
    return;
  }

  void CmdParser::addInputs(const std::string & param)
  {
    YAML::Node rootNode;
    if (!rootNode)
    {
      std::cerr << "Root node is empty\n";
      exit(EXIT_FAILURE);
    }

    rootNode["inputs"][param];
    return;
  }

  void CmdParser::addOutputs(const std::string & param)
  {
    YAML::Node rootNode;
    if (!rootNode)
    {
      std::cerr << "Root node is empty\n";
      exit(EXIT_FAILURE);
    }

    rootNode["outputs"][param];
    return;
  }

  std::string CmdParser::checkDefault(const std::string & param)
  {
    YAML::Node input;
    std::string defaultValue = "";
    if (input[param]["default"]) {
      defaultValue = (input[param]["default"]).as<std::string>();;
    }
    return defaultValue;
  }

  void CmdParser::cwlrunner(const std::string & cwl_spec_path, const std::string & cwl_input_path, bool getDefaultFlag)
  {
    std::string cmd = "cwl-runner ";
    if (!getDefaultFlag) {
      cmd += cwl_spec_path + " " + cwl_input_path;
    }
    cmd += cwl_spec_path;
    std::system(cmd.c_str());
    //ShellExecute(0, "open", "cmd.exe", cmd.c_str(), 0, SW_HIDE);;
  }

  std::string CmdParser::getLaconic(const std::string & execParamToCheck)
  {
    int noMoreChecks = 0; // ensures that extra checks are not done for parameters
    if (execParamToCheck.empty())
    {
      std::cerr << "Parameter cannot be an empty string. Please try again.\n";
      exit(EXIT_FAILURE);
    }
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    size_t i = 0;
    while ((i < m_requiredParameters.size()) && (noMoreChecks < 1))
    {
      if ((m_requiredParameters[i].verbose == execParamToCheck))
      {
        return m_requiredParameters[i].laconic;
        noMoreChecks = 1;
      }
      i++;
    }

    i = 0;
    while ((i < m_optionalParameters.size()) && (noMoreChecks < 1))
    {
      if ((m_optionalParameters[i].verbose == execParamToCheck))
      {
        return m_optionalParameters[i].laconic;
        noMoreChecks = 1;
      }
      i++;
    }
    return "";
  }

  void CmdParser::logCWL(const std::string &inpFileName, const std::string &cwlFileName) 
  {  
    ////create a timecode specific folder
    //// current date/time based on current system
    //time_t now = time(0);

    //struct tm timeinfo;
    //localtime_s(&timeinfo, &now);

    //// print various components of tm structure.
    //std::cout << "Year: " << 1900 + timeinfo.tm_year << std::endl;
    //int Y = 1900 + timeinfo.tm_year;
    //std::cout << "Month: " << 1 + timeinfo.tm_mon << std::endl;
    //int M = 1 + timeinfo.tm_mon;
    //std::cout << "Day: " << timeinfo.tm_mday << std::endl;
    //int D = timeinfo.tm_mday;
    //std::cout << "Time: " << 1 + timeinfo.tm_hour << ":";
    //int H = timeinfo.tm_hour;
    //std::cout << 1 + timeinfo.tm_min << ":";
    //int Mi = timeinfo.tm_min;
    //std::cout << 1 + timeinfo.tm_sec << std::endl;
    //int S = 1 + timeinfo.tm_sec;

    //std::string timestampStr;
    //std::stringstream convD, convM, convY, convH, convMi, convS;
    //convD << D;
    //convM << M;
    //convY << Y;
    //convH << H;
    //convMi << Mi;
    //convS << S;
    //std::cout << "Timestamp:" << std::endl;
    //timestampStr = convD.str() + '.' + convM.str() + '.' + convY.str() + '-' + convH.str() + ':' + convMi.str() + ':' + convS.str();
    //std::cout << timestampStr << std::endl;

    //namespace fs = std::experimental::filesystem;
    //fs::path inpSourceFile = inpFileName;
    //fs::path cwlSourceFile = cwlFileName;
    //fs::path targetParent = "m_exeName/" + timestampStr;
    //auto target1 = targetParent / inpSourceFile.filename(); // sourceFile.filename() returns "sourceFile.ext".
    //auto target2 = targetParent / cwlSourceFile.filename();
    //try // If you want to avoid exception handling, then use the error code overload of the following functions.
    //{
    //  fs::create_directories(targetParent); // Recursively create target directory if not existing.
    //  fs::copy_file(inpSourceFile, target1, fs::copy_options::overwrite_existing);
    //  fs::create_directories(targetParent); // Recursively create target directory if not existing.
    //  fs::copy_file(cwlSourceFile, target2, fs::copy_options::overwrite_existing);
    //}
    //catch (std::exception& e) // Not using fs::filesystem_error since std::bad_alloc can throw too.  
    //{
    //  std::cout << e.what();
    //}    
  }


  void CmdParser::writeConfigFile(const std::string &dirName)
  {
    if (!checkMaxLen)
    {
      getMaxLength();
    }

    std::string dirName_wrap;
    if (!cbica::directoryExists(dirName) || (dirName == ""))
    {
      dirName_wrap = cbica::makeTempDir();
    }
    dirName_wrap = cbica::stringReplace(dirName, "\\", "/");
    if (dirName_wrap.substr(dirName_wrap.length() - 1) != "/")
    {
      dirName_wrap += "/";
    }

    std::string fileName = dirName_wrap + m_exeName + ".txt";



    //#if (_WIN32)
    //    if (_access(fileName.c_str(), 6) == -1)
    //    {
    //      std::cerr << "No write permission for the specified config file.\n";
    //      exit(EXIT_FAILURE);
    //    }
    //#else
    //    if (access(fileName.c_str(), R_OK && W_OK) != 0)
    //    {
    //      std::cerr << "No write permission for the specified config file.\n";
    //      exit(EXIT_FAILURE);
    //    }
    //#endif

    std::ofstream file;
    file.open(fileName.c_str());

    if (file.is_open())
    {
      for (size_t i = 0; i < m_requiredParameters.size(); i++)
      {
        file << getSeparator(Param) << m_requiredParameters[i].verbose << getSeparator(Param) <<
          " " << getSeparator(DataType) << m_requiredParameters[i].dataType_string << getSeparator(DataType) <<
          " " << getSeparator(DataRange) << m_requiredParameters[i].dataRange << getSeparator(DataRange) <<
          " " << m_requiredParameters[i].descriptionLine1 + " " + m_requiredParameters[i].descriptionLine2 + " " +
          m_requiredParameters[i].descriptionLine3 + " " + m_requiredParameters[i].descriptionLine4 + " " +
          m_requiredParameters[i].descriptionLine5 << "\n";
      }

      for (size_t i = 0; i < m_optionalParameters.size(); i++)
      {
        file << getSeparator(Param) << m_optionalParameters[i].verbose << getSeparator(Param) <<
          " " << getSeparator(DataType) << m_optionalParameters[i].dataType_string << getSeparator(DataType) <<
          " " << getSeparator(DataRange) << m_optionalParameters[i].dataRange << getSeparator(DataRange) <<
          " " << m_optionalParameters[i].descriptionLine1 + " " + m_optionalParameters[i].descriptionLine2 + " " +
          m_optionalParameters[i].descriptionLine3 + " " + m_optionalParameters[i].descriptionLine4 + " " +
          m_optionalParameters[i].descriptionLine5 << "\n";
      }
    }

    file.close();

    //std::cout << "Config file written with path: '" << fileName << "'\n";
    return;
  }

  std::vector< Parameter > CmdParser::readConfigFile(const std::string &path_to_config_file, bool getDescription)
  {
    std::vector< Parameter > returnVector;
    std::ifstream inputFile(path_to_config_file.c_str());
    if (!inputFile)
    {
      std::cerr << "File '" << path_to_config_file << "' not found.\n";
      exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(inputFile, line))
    {
      std::string parameter, parameterDataType, parameterDataRange, parameterDescription = "";
      for (size_t i = 0; i < line.length(); i++)
      {
        parameter = cbica::iterateOverStringAndSeparators(line, i,
#ifdef _WIN32
          Separator::
#endif
          Param);
        i = i + 2;
        parameterDataType = cbica::iterateOverStringAndSeparators(line, i,
#ifdef _WIN32
          Separator::
#endif
          DataType);
        i = i + 2;
        parameterDataRange = cbica::iterateOverStringAndSeparators(line, i,
#ifdef _WIN32
          Separator::
#endif
          DataRange);
        if (getDescription)
        {
          i = i + 1;
          parameterDescription = cbica::iterateOverStringAndSeparators(line, i, 10);
        }
        i = line.length();
        returnVector.push_back(Parameter("", parameter, parameterDataType, parameterDataRange, parameterDescription));
      }
    }

    inputFile.close();
    return returnVector;
  }

}
