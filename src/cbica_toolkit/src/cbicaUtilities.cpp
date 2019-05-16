/**
\file  cbicaUtilities.cpp

\brief Some basic utility functions.

Dependecies: OpenMP

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#if (_WIN32)
#define NOMINMAX
#include <direct.h>
#include <windows.h>
#include <conio.h>
#include <lmcons.h>
#include <Shlobj.h>
#include <filesystem>
#include <psapi.h>
#define GetCurrentDir _getcwd
bool WindowsDetected = true;
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
#if (__APPLE__)
  #include <mach-o/dyld.h>
  #include <sys/sysctl.h>
  #include <mach/vm_statistics.h>
  #include <mach/mach_types.h>
  #include <mach/mach_init.h>
  #include <mach/mach_host.h>
#else
  #include <sys/sysinfo.h>
#endif
#define GetCurrentDir getcwd
bool WindowsDetected = false;
static const char  cSeparator = '/';
//  static const char* cSeparators = "/";
#endif

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
#ifdef _OPENMP
#include <omp.h>
#endif
#include <thread>

#include "cbicaUtilities.h"
#include "yaml-cpp/yaml.h"

namespace cbica
{
  //====================================== Folder stuff ====================================//

  bool fileExists(const std::string &fName)
  {
    std::ifstream file_exists(fName.c_str());
    if (file_exists.good())
      return true;
    else
      return false;
  }

  bool directoryExists(const std::string &dName)
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

  bool isFile(const std::string &path)
  {
    return cbica::fileExists(path);
  }

  bool isDir(const std::string &path)
  {
    return cbica::directoryExists(path);
  }

  bool exists(const std::string &path)
  {
    struct stat info;

    if (stat(path.c_str(), &info) != 0)
      return false;
    else if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows
      return true;
    else
      return true;
  }

  std::vector<std::string> getCWLFilesInApplicationDir() 
  {
    auto appDir = getExecutablePath();

#ifdef __APPLE__
    appDir += "../Resources/bin";
#endif

    auto filesInDir = filesInDirectory(appDir);
    auto cwlFiles = filesInDir;
    cwlFiles.clear();
    for (size_t i = 0; i < filesInDir.size(); i++)
    {
      if (getFilenameExtension(filesInDir[i], false).find(".cwl") != std::string::npos)
      {
        cwlFiles.push_back(filesInDir[i]);
      }
    }
    //std::vector<std::string> files;

    //#ifdef _WIN32
    //  WIN32_FIND_DATA data;
    //  HANDLE hFind = FindFirstFile("\\*", &data);

    //  if ( hFind != INVALID_HANDLE_VALUE ) {
    //    do {
    //      files.push_back(data.cFileName);
    //    } while (FindNextFile(hFind, &data));
    //    FindClose(hFind);
    //  }
    //#else
    //  DIR *dir;
    //  struct dirent *ent;
    //  if ((dir = opendir(".")) != NULL) {
    //    /* print all the files and directories within directory */
    //    while ((ent = readdir (dir)) != NULL) {
    //      if (ent->d_type == DT_REG) {  
    //        files.push_back(ent->d_name);
    //      }
    //    }
    //    closedir (dir);
    //  } else {
    //    /* could not open directory */
    //    perror ("");
    //    return files;
    //  }
    //#endif

    //// Prune non cwl files
    //std::vector<std::string> cwlfiles;
    //for(auto const& value: files) {

    //  if (value.substr(value.size() - 4) == ".cwl") {
    //    cwlfiles.push_back(value);
    //  }

    //}

    //// Sort cwl files
    //std::sort(cwlfiles.begin(), cwlfiles.end());

    return cwlFiles;

  }

  std::string getEnvironmentVariableValue(const std::string &environmentVariable)
  {
    std::string returnString = "";
    char tempValue[FILENAME_MAX];
#if defined(_WIN32)
    char tmp[FILENAME_MAX];
    size_t size = FILENAME_MAX;
    getenv_s(&size, tmp, size, environmentVariable.c_str()); // does not work, for some reason - needs to be tested
    std::string temp = cbica::replaceString(tmp, "\\", "/");
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

  std::string getCurrentProcessID()
  {
    return std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()));
  }

  std::string createTmpDir()
  {
    std::string returnDir = "", tempCheck, homeEnv;
#if defined(_WIN32)
    homeEnv = "USERPROFILE";
#else
    homeEnv = "HOME";
#endif

    auto exeTempDir = cbica::getEnvironmentVariableValue(homeEnv) + "/.cbicaTemp/"/* + cbica::getExecutableName()*/;
    if (!directoryExists(exeTempDir))
    {
      createDir(exeTempDir);
      auto temp = cbica::stringSplit(cbica::getCurrentLocalTime(), ":");
    }

    tempCheck = exeTempDir + "/tmp_" + getCurrentProcessID();

    if (cbica::directoryExists(tempCheck))
    {
      auto temp = cbica::stringSplit(cbica::getCurrentLocalTime(), ":");
      tempCheck += temp[0] + temp[1] + temp[2] + "/";
    }

    if (isDir(tempCheck))
    {
      for (size_t i = 1; i <= FILENAME_MAX; i++)
      {
        returnDir = tempCheck + std::to_string(i);
        if (!isDir(returnDir))
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

  std::string createTemporaryDirectory()
  {
    return createTmpDir();
  }

  std::string makeTemporaryDirectory()
  {
    return createTmpDir();
  }

  std::string makeTempDir()
  {
    return createTmpDir();
  }

  bool createDir(const std::string &dir_name)
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

  bool makeDir(const std::string &dir_name)
  {
    return createDir(dir_name);
  }

  bool createDirectory(const std::string &dir_name)
  {
    return createDir(dir_name);
  }

  bool makeDirectory(const std::string &dir_name)
  {
    return createDir(dir_name);
  }

  bool createFolder(const std::string &dir_name)
  {
    return createDir(dir_name);
  }

  int removeDirectoryRecursively(const std::string &dirname, bool bDeleteSubdirectories = true)
  {
    if (!directoryExists(dirname))
    {
      std::cerr << "Supplied directory name wasn't found: " << dirname << std::endl;
      exit(EXIT_FAILURE);
    }
#if defined(_WIN32)
    bool bSubdirectory = false;       // Flag, indicating whether
    // subdirectories have been found
    HANDLE hFile;                     // Handle to directory
    std::string strFilePath;          // Filepath
    std::string strPattern;           // Pattern
    WIN32_FIND_DATA FileInformation;  // File information

    strPattern = dirname + "/*.*";
    hFile = ::FindFirstFile(strPattern.c_str(), &FileInformation);
    if (hFile != INVALID_HANDLE_VALUE)
    {
      do
      {
        if (FileInformation.cFileName[0] != '.')
        {
          strFilePath.erase();
          strFilePath = dirname + "/" + FileInformation.cFileName;

          if (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
          {
            if (bDeleteSubdirectories)
            {
              // Delete subdirectory
              int iRC = cbica::removeDirectoryRecursively(strFilePath, bDeleteSubdirectories);
              if (iRC)
                return iRC;
            }
            else
              bSubdirectory = true;
          }
          else
          {
            // Set file attributes
            if (::SetFileAttributes(strFilePath.c_str(),
              FILE_ATTRIBUTE_NORMAL) == FALSE)
              return ::GetLastError();

            // Delete file
            if (::DeleteFile(strFilePath.c_str()) == FALSE)
              return ::GetLastError();
          }
        }
      } while (::FindNextFile(hFile, &FileInformation) == TRUE);

      // Close handle
      ::FindClose(hFile);

      DWORD dwError = ::GetLastError();
      if (dwError != ERROR_NO_MORE_FILES)
        return dwError;
      else
      {
        if (!bSubdirectory)
        {
          // Set directory attributes
          if (::SetFileAttributes(dirname.c_str(),
            FILE_ATTRIBUTE_NORMAL) == FALSE)
            return ::GetLastError();

          // Delete directory
          if (::RemoveDirectory(dirname.c_str()) == FALSE)
            return ::GetLastError();
        }
      }
    }

    return 0;
#else
    std::string passString = "rm -rf " + dirname;
    if (std::system(passString.c_str()) != 0)
      std::cerr << "Error during delete.\n";
#endif
    return 0;
  }

  bool removeDir(const std::string &path)
  {
    return deleteDir(path);
  }

  bool deleteDir(const std::string &path)
  {
    //if (removeDirectoryRecursively(path, true) != 0)
    //  return false;
    //return true;
    //return removeDirectoryRecursively(path);
#if defined(_WIN32)
    if (_rmdir(path.c_str()) == -1)
      return false;
#else
    std::string passString = "rmdir " + path;
    if (system(passString.c_str()) != 0)
      return false;
    #endif

    return true;
  }

  bool copyDir(const std::string &inputFolder, const std::string &destination, bool recursion)
  {
    if (cbica::isDir(inputFolder))
    {
      if (!cbica::isDir(destination))
      {
        cbica::createDir(destination);
      }
#ifdef _WIN32
      WIN32_FIND_DATA FindFileData;
      HANDLE hFind;
      const size_t MAX_FILE_SIZE = 1025;

      std::string source_tmp = inputFolder;
      // ensure that all folders end with "/"
      if (source_tmp[source_tmp.length() - 1] != '/')
      {
        source_tmp.append("/");
      }

      // check for all files inside
      source_tmp.append("*");

      hFind = FindFirstFile(source_tmp.c_str(), &FindFileData);
      if (hFind == NULL)
      {
        return false;
      }
      do
      {
        if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
        {
          if (strcmp(FindFileData.cFileName, "."))
          {
            if (strcmp(FindFileData.cFileName, ".."))
            {
              std::string newSource = inputFolder + std::string(FindFileData.cFileName) + "/",
                newDest = destination + std::string(FindFileData.cFileName) + "/";
              createDir(newDest);
              //CreateDirectory(newDest.c_str(), NULL);
              copyDir(newSource, newDest);
            }
          }
        }
        else
        {
          std::string sourceFile = inputFolder + std::string(FindFileData.cFileName), destinationFile = destination + std::string(FindFileData.cFileName);
          copyFile(sourceFile, destinationFile);
          //BOOL l_bRet = CopyFile(sourceFile.c_str(), destinationFile.c_str(), TRUE);
        }

      } while (FindNextFile(hFind, &FindFileData));
      FindClose(hFind);

#else
      //if (recursion)
      //{
      //  recur = "-r ";
      //}
      //system(std::string("cp " + recur + " " + inputFolder + " " + destination).c_str());

      DIR *dir = opendir(inputFolder.c_str());                //Assuming absolute pathname here.
      if (dir)
      {
        char Path[256], *EndPtr = Path;
        struct dirent *e;
        strcpy(Path, inputFolder.c_str());                  //Copies the current path to the 'Path' variable.
        EndPtr += strlen(inputFolder.c_str());              //Moves the EndPtr to the ending position.
        while ((e = readdir(dir)) != NULL) //Iterates through the entire directory
        {
          struct stat info;                //Helps us know about stuff
          strcpy(EndPtr, e->d_name);       //Copies the current filename to the end of the path, overwriting it with each loop.
          if (!stat(Path, &info)) //stat returns zero on success.
          {
            if (S_ISDIR(info.st_mode)) //Are we dealing with a directory?
            {
              //Make corresponding directory in the target folder here.
              copyDir(Path, EndPtr);   //Calls this function AGAIN, this time with the sub-name.
            }
            else if (S_ISREG(info.st_mode)) //Or did we find a regular file?
            {
              copyFile(Path, EndPtr);
            }
          }
        }
      }
#endif

      if (cbica::isDir(destination))
      {
        return true;
      }
      else
      {
        std::cerr << "Something went wrong when trying to do the copy. Ensure you have write access to destination.\n";
        return false;
      }
    }
    else
    {
      std::cerr << "The input folder '" << inputFolder << "' cannot be verified as a folder.\n";
      return false;
    }
  }

  bool copyDirectory(const std::string &inputFolder, const std::string &destination, bool recursion)
  {
    return copyDir(inputFolder, destination, recursion);
  }

  bool copyFolder(const std::string &inputFolder, const std::string &destination, bool recursion)
  {
    return copyDir(inputFolder, destination, recursion);
  }

  bool copyFile(const std::string &inputFile, const std::string &destination)
  {
    if (cbica::fileExists(inputFile))
    {
      std::ifstream src(inputFile.c_str(), std::ios::binary);
      std::ofstream dst(destination.c_str(), std::ios::binary);

      dst << src.rdbuf();

      if (cbica::fileExists(destination))
      {
        return true;
      }
      else
      {
        std::cerr << "Something went wrong when trying to do the copy. Ensure you have write access to destination.\n";
        return false;
      }
    }
    else
    {
      std::cerr << "The input file '" << inputFile << "' cannot be verified as a file.\n";
      return false;
    }
  }

  size_t getFileSize(const std::string &inputFile)
  {
    std::streampos begin, end;
    std::ifstream myfile(inputFile.c_str(), std::ios::binary);
    begin = myfile.tellg();
    myfile.seekg(0, std::ios::end);
    end = myfile.tellg();
    myfile.close();

    return (end - begin);

    /* // if using filesystem
    std::tr2::sys::path folderPath(inputFile);
    return file_size(filePath);
    */
  }

  bool IsCompatible(const std::string inputVersionFile)
  {
    auto config = YAML::LoadFile(inputVersionFile);

    auto currentCollectionVersion = std::stoi(cbica::replaceString(config["Version"].as< std::string >().c_str(), ".", "").c_str());
    auto minimumVersion = std::stoi(cbica::replaceString(config["Minimum"].as< std::string >().c_str(), ".", "").c_str());
    auto maximumVersion = std::stoi(cbica::replaceString(config["Maximum"].as< std::string >().c_str(), ".", "").c_str());
    auto currentPackageVersion = std::stoi(cbica::replaceString(std::string(PROJECT_VERSION), ".", "").c_str());

    if (currentPackageVersion == currentCollectionVersion)
    {
      return true;
    }
    if (currentPackageVersion < minimumVersion)
    {
      return false;
    }
    if (currentPackageVersion > maximumVersion)
    {
      return false;
    }

    return true;
  }

  size_t getFolderSize(const std::string &rootFolder)
  {
    size_t f_size = 0;
#if _WIN32
    std::experimental::filesystem::path folderPath(rootFolder);
    if (exists(folderPath))
    {
		std::experimental::filesystem::directory_iterator end_itr;
      for (std::experimental::filesystem::directory_iterator dirIte(rootFolder); dirIte != end_itr; ++dirIte)
      {
		  std::experimental::filesystem::path filePath;
#if (_MSC_VER >= 1900)
        filePath = std::experimental::filesystem::system_complete(dirIte->path());
#else
        filePath = complete(dirIte->path(), folderPath);
#endif
        if (!is_directory(dirIte->status()))
        {
          f_size += file_size(filePath);
        }
        else
        {
          f_size += getFolderSize(filePath.string());
        }
      }
    }
    else
    {
      std::cerr << "Folder not found.\n";
      exit(EXIT_FAILURE);
    }
#else
    DIR *d;
    struct dirent *de;
    struct stat buf;
    int exists;

    d = opendir(".");
    if (d == NULL)
    {
      perror("prsize");
      exit(1);
    }

    f_size = 0;

    for (de = readdir(d); de != NULL; de = readdir(d))
    {
      exists = stat(de->d_name, &buf);
      if (exists < 0)
      {
        fprintf(stderr, "Couldn't stat %s\n", de->d_name);
      }
      else
      {
        f_size += buf.st_size;
      }
    }
    closedir(d);
#endif
    return f_size;
  }

  size_t getDirSize(const std::string &rootFolder)
  {
    return getFolderSize(rootFolder);
  }

  size_t getDirectorySize(const std::string &rootFolder)
  {
    return getFolderSize(rootFolder);
  }

  //======================================== OS stuff ======================================//

  std::string getFilenameBase(const std::string &filename, bool checkFile)
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

  std::string getFilenameExtension(const std::string &filename, bool checkFile)
  {
    if (checkFile)
    {
      if (!fileExists(filename))
      {
        std::cerr << "[getFilenameExtension()] Supplied file name'" << filename << "'wasn't found.\n";
        exit(EXIT_FAILURE);
      }
    }
    std::string path, base, ext;
    splitFileName(filename, path, base, ext);

    return ext;
  }

  std::string getFilenamePath(const std::string &filename, bool checkFile)
  {
    if (directoryExists(filename))
    {
      return filename;
    }
    else
    {
      if (checkFile)
      {
        if (!fileExists(filename))
        {
          std::cerr << "[getFilenamePath()] Supplied file name'" << filename << "'wasn't found.\n";
          exit(EXIT_FAILURE);
        }
      }
      std::string path, base, ext;
      splitFileName(filename, path, base, ext);

      return path;
    }
  }

  std::string getExecutableName()
  {
    std::string return_string;
#if defined(_WIN32)
    //! Initialize pointers to file and user names
    char filename[FILENAME_MAX];
    GetModuleFileNameA(NULL, filename, FILENAME_MAX);
    std::string path, ext;
    splitFileName(filename, path, return_string, ext);
    filename[0] = '\0';
    //_splitpath_s(filename, NULL, NULL, NULL, NULL, filename, NULL, NULL, NULL);
#elif __APPLE__
    char path[PATH_MAX];
    uint32_t size = PATH_MAX - 1;
    if (path != NULL)
    {
      if (_NSGetExecutablePath(path, &size) == 0)
      {
        return_string = getFilenameBase(std::string(path));
      }
    }
#else
    return_string = getEnvironmentVariableValue("_");
    return_string = cbica::replaceString(return_string, "./", "");
    char path[PATH_MAX];
    if (path != NULL)
    {
      auto ret = ::readlink("/proc/self/exe", path, sizeof(path)-1);
      if (ret == -1)
      {
        //free(path);
        path[0] = '\0';
      }
      path[ret] = 0;
    }
    auto temp = std::string(path);
    return_string = getFilenameBase(temp);
    path[0] = '\0';
#endif

    return return_string;
  }

  std::string getExecutablePath()
  {
    std::string path, base, ext;
    splitFileName(getFullPath(), path, base, ext);
  // #ifdef __APPLE__
  //   return "/Applications/CaPTk_1.6.2.Beta.app/Contents/MacOS/";
  // #endif
    return cbica::normPath(path) + "/";
  }

  std::string getFullPath()
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

  std::string getUserName()
  {
#if defined(_WIN32)
    //! Initialize pointers to file and user names
    char username[FILENAME_MAX];
    DWORD username_len = FILENAME_MAX;
    GetUserName(username, &username_len);
#else
    //! Initialize pointers to file and user names
    char username[FILENAME_MAX];
    size_t username_len = FILENAME_MAX;
    getlogin_r(username, username_len);
#endif

    std::string return_string = std::string(username);
    username[0] = '\0';

    return return_string;
  }

  std::string getUserHomeDirectory()
  {
    std::string homeProfile;
#if defined(_WIN32)
    homeProfile = "USERPROFILE";
#else
    homeProfile = "HOME";
#endif

    return cbica::getEnvironmentVariableValue(homeProfile);
  }

  std::string getCWD()
  {
    std::string wd = "";

    // use c++ convention
    char* buffer = GetCurrentDir(NULL, 0);

    if (buffer)
    {
      wd = buffer;
      buffer[0] = '\0';
    }
    wd = replaceString(wd, "\\", "/");
    if (wd[wd.length() - 1] != '/')
    {
      wd.append("/");
    }

    return wd;
  }

  //! Check for separators [Internal_Function]
  inline bool issep(char c)
  {
#if defined(_WIN32)
    return c == '/' || c == '\\';
#else
    return c == '/';
#endif
  }

  //! Check for absolute path [Internal_Function]
  bool isabs(const std::string& path)
  {
    size_t i = 0;
#if defined(_WIN32)
    if (path.size() > 1 && path[1] == ':') i = 2;
#endif
    return i < path.size() && issep(path[i]);
  }

  //! Joins separators based on OS [Internal_Function]
  std::string join(const std::string& base, const std::string& path)
  {
    if (base.empty() || isabs(path))
      return path;
    if (issep(base[base.size() - 1]))
      return (base + path);

#if defined(_WIN32)
    return base + '\\' + path;
#else
    return base + '/' + path;
#endif
  }

  std::string normPath(const std::string &path)
  {
    if (path.empty())
      return "";
    char drive[3] = { '\0', ':', '\0' };
    size_t i = 0;
#if defined(_WIN32)
    if (path.size() > 1 && path[1] == ':')
    {
      drive[0] = path[0];
      i = 2;
    }
#endif
    std::string norm_path = drive;
    bool abs = issep(path[i]);
    if (abs)
    {
#if defined(_WIN32)
      while (i <= path.size() && issep(path[i]))
      {
        norm_path += cSeparator;
        i++;
      }
#else
      norm_path += cSeparator;
#endif
    }
    std::string current;
    std::vector<std::string> parts;
    while (i <= path.size())
    {
      if (issep(path[i]) || path[i] == '\0')
      {
        if (current == "..")
        {
          if (!abs && (parts.empty() || parts.back() == ".."))
          {
            parts.push_back(current);
          }
          else if (!parts.empty())
          {
            parts.pop_back();
          }
        }
        else if (current != "" && current != ".")
        {
          parts.push_back(current);
        }
        current.clear();
      }
      else
      {
        current += path[i];
      }
      i++;
    }
    for (i = 0; i < parts.size(); i++)
    {
      norm_path = join(norm_path, parts[i]);
    }
    std::replace(norm_path.begin(), norm_path.end(), '\\', '/');
    return norm_path.empty() ? "." : norm_path;
  }

  std::string normalizePath(const std::string &path)
  {
    return normPath(path);
  }

  //! Absolute path [Internal_Function]
  std::string absPath(const std::string &path)
  {
    return normPath(join(getCWD(), path));
  }

  //! Split Drive name from path [Internal_Function]
  inline void splitDrive(const std::string& path, std::string& drive, std::string& tail)
  {
#if defined(_WIN32)
    if (path.size() > 1 && path[1] == ':')
    {
      tail = path.substr(2);
      drive = path[0]; drive += ':';
    }
    else
#endif
    {
      tail = path;
      drive = "";
    }
  }

  //! Split Drive name from path [Internal_Function]
  inline std::vector<std::string> splitDrive(const std::string& path)
  {
    std::vector<std::string> parts(2, "");
    splitDrive(path, parts[0], parts[1]);
    return parts;
  }

  std::string relPath(const std::string &path, const std::string &base)
  {
    // if relative path is given just return it
    if (!isabs(path))
      return path;
    // normalize paths
    std::string norm_path = normPath(path);
    std::string norm_base = normPath(join(getCWD(), base));
    // check if paths are on same drive
#if defined(_WIN32)
    std::string drive = splitDrive(norm_path)[0];
    std::string base_drive = splitDrive(norm_base)[0];
    if (drive != base_drive)
    {
      std::cerr << "Error: Path is on drive " << drive << ", start is on drive " << base_drive;
    }
#endif
    // find start of first path component in which paths differ
    std::string::const_iterator b = norm_base.begin();
    std::string::const_iterator p = norm_path.begin();
    size_t pos = 0;
    size_t i = 0;
    while (b != norm_base.end() && p != norm_path.end())
    {
      if (issep(*p))
      {
        if (!issep(*b))
          break;
        pos = i;
      }
      else if (*b != *p)
      {
        break;
      }
      b++; p++; i++;
    }
    // set pos to i (in this case, the size of one of the paths) if the end
    // of one path was reached, but the other path has a path separator
    // at this position, this is required below
    if ((b != norm_base.end() && issep(*b)) || (p != norm_path.end() && issep(*p)))
      pos = i;
    // skip trailing separator of other path if end of one path reached
    if (b == norm_base.end() && p != norm_path.end() && issep(*p))
      p++;
    if (p == norm_path.end() && b != norm_base.end() && issep(*b))
      b++;
    // if paths are the same, just return a period (.)
    //
    // Thanks to the previous skipping of trailing separators, this condition
    // handles all of the following cases:
    //
    //    base := "/usr/bin"  path := "/usr/bin"
    //    base := "/usr/bin/" path := "/usr/bin/"
    //    base := "/usr/bin"  path := "/usr/bin/"
    //    base := "/usr/bin/" path := "/usr/bin"
    if (b == norm_base.end() && p == norm_path.end())
      return ".";
    // otherwise, pos is the index of the last slash for which both paths
    // were identical; hence, everything that comes after in the original
    // path is preserved and for each following component in the base path
    // a "../" is prepended to the relative path
    std::string rel_path;
    // truncate base path with a separator as for each "*/" path component,
    // a "../" will be prepended to the relative path
    if (b != norm_base.end() && !issep(norm_base[norm_base.size() - 1]))
    {
      // attention: This operation may invalidate the iterator b!
      //            Therefore, remember position of iterator and get a new one.
      size_t pos = b - norm_base.begin();
      norm_base += cSeparator;
      b = norm_base.begin() + pos;
    }
    while (b != norm_base.end())
    {
      if (issep(*b))
      {
        rel_path += "..";
        rel_path += cSeparator;
      }
      b++;
    }
    if (pos + 1 < norm_path.size())
      rel_path += norm_path.substr(pos + 1);
    // remove trailing path separator
    if (issep(rel_path[rel_path.size() - 1]))
    {
      rel_path.erase(rel_path.size() - 1);
    }
    return rel_path;
  }

  std::string relativePath(const std::string &path, const std::string &base)
  {
    return relPath(path, base);
  }

  std::string realPath(const std::string &path)
  {
    std::string curr_path = join(getCWD(), path);
#if defined(_WIN32)
    // nothing extra required
#else
    char *actualPath = realpath(const_cast<char *>(curr_path.c_str()), NULL);
    curr_path = std::string(actualPath);
    /*
    // use stringstream and std::getline() to split absolute path at slashes (/)
    std::stringstream ss(curr_path);
    curr_path.clear();
    std::string fname;
    std::string prev_path;
    std::string next_path;
    char slash;
    ss >> slash; // root slash
    while( getline(ss, fname, '/') )
    {
    // current absolute path
    curr_path += '/';
    curr_path += fname;
    // if current path is a symbolic link, follow it
    if( isLink(curr_path) )
    {
    // for safety reasons, restrict the depth of symbolic links followed
    for( unsigned int i = 0; i < 100; i++ )
    {
    char *buffer=NULL, *newbuf=NULL;
    size_t buflen = 256;
    for(;;)
    {
    newbuf = reinterpret_cast<char*>(realloc(buffer, buflen * sizeof(char)) );
    if( !newbuf )
    break;
    buffer = newbuf;
    int n = ::newlink(path.c_str(), buffer, buflen);
    if( n<0 )
    break;
    if( static_cast<size_t>(n)<buflen )
    {
    buffer[n] = '\0';
    next_path = buffer;
    break;
    }
    buflen+=256;
    }
    free(buffer);
    //next_path = os::readlink(curr_path);
    if( next_path.empty() )
    {
    // if real path could not be determined because of permissions
    // or invalid path, return the original path
    break;
    }
    else
    {
    curr_path = join(prev_path, next_path);
    if( !isLink(next_path) )
    break;
    }
    }
    // if real path could not be determined with the given maximum number
    // of loop iterations (endless cycle?) or one of the symbolic links
    // could not be read, just return original path as absolute path
    if( isLink(next_path) )
    return absPath(path);
    }
    // memorize previous path used as base for abspath()
    prev_path = curr_path;
    }
    */
#endif
    // normalize path after all symbolic links were resolved
    return normPath(curr_path);
  }

  bool isLink(const std::string &path)
  {
#if defined(_WIN32)
    std::cerr << "Windows doesn't support ways to distinguish between hard and soft links.\n";
    return false;
#else
    struct stat info;
    if (lstat(path.c_str(), &info) != 0)
      return false;
    return S_ISLNK(info.st_mode);
#endif
  }

  bool isSymbolicLink(const std::string &path)
  {
    return isLink(path);
  }

  bool makeSymbolicLink(const std::string &input_fileName, const std::string &ouput_fileName)
  {
    if (!cbica::fileExists(input_fileName))
    {
      std::cerr << "Supplied file name wasn't found.\n";
      exit(EXIT_FAILURE);
    }
#if defined(_WIN32)
    if (IsUserAnAdmin())
    {
      if (CreateSymbolicLink(input_fileName.c_str(), ouput_fileName.c_str(), 0) != 0)
        return true;
      else
        return false;
    }
    else
    {
      std::cerr << "Windows doesn't let non-admins create soft links.\n";
      return false;
    }
#else
    if (symlink(input_fileName.c_str(), ouput_fileName.c_str()) == 0)
      return true;
    else
      return false;
#endif
  }

  bool setEnvironmentVariable(const std::string &variable_name, const std::string &variable_value)
  {
    std::string totalVariable = variable_name + "=" + variable_value;
    try
    {
#if defined(_WIN32)
      int test = _putenv(totalVariable.c_str());
#else
      putenv(cbica::constCharToChar(totalVariable));
#endif
    }
    catch (const std::exception &e)
    {
      std::cerr << "Exception caught: " << e.what() << std::endl;
      return false;
    }

    return true;
  }

  bool deleteEnvironmentVariable(const std::string &variable_name)
  {
    return cbica::setEnvironmentVariable(variable_name, "");
  }

  std::vector< std::string > filesInDirectory(const std::string &dirName, bool returnFullPath)
  {
    if (!cbica::directoryExists(dirName))
    {
      std::cerr << "Supplied directory name wasn't found: " << dirName << std::endl;
      exit(EXIT_FAILURE);
    }
    std::vector< std::string > allFiles;
    std::string dirName_wrap = cbica::normPath(dirName);
    if (dirName_wrap[dirName_wrap.length() - 1] != '/')
    {
      dirName_wrap.append("/");
    }
#if defined(_WIN32)
    {
      char* search_path = cbica::constCharToChar((dirName_wrap + "*.*").c_str());
      WIN32_FIND_DATA fd;
      HANDLE hFind = ::FindFirstFile(search_path, &fd);
      if (hFind != INVALID_HANDLE_VALUE)
      {
        do
        {
          if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
          {
            if (returnFullPath)
            {
              allFiles.push_back(dirName_wrap + fd.cFileName);
            }
            else
            {
              allFiles.push_back(fd.cFileName);
            }
          }
        } while (::FindNextFile(hFind, &fd));
        ::FindClose(hFind);
      }
    }
#else
    {
      DIR *dp;
      struct dirent *dirp;
      if ((dp = opendir(dirName.c_str())) == NULL)
      {
        std::cerr << "Error(" << errno << ") occurred while opening directory '" <<
          dirName << "'\n";
      }

      while ((dirp = readdir(dp)) != NULL)
      {
        auto returnBase = std::string(dirp->d_name);
        if((returnBase.compare(".") != 0) && (returnBase.compare("..") != 0) && (returnBase.compare("~") != 0))
        {
          if (returnFullPath)
          {
            allFiles.push_back(dirName_wrap + std::string(dirp->d_name));
          }
          else
          {
            allFiles.push_back(std::string(dirp->d_name));
          }
        }
      }
      closedir(dp);
    }
#endif

  std::sort(allFiles.begin(), allFiles.end());
  return allFiles;
  }

  std::vector<std::string> subdirectoriesInDirectory(const std::string &dirName, bool recursiveSearch, bool returnFullPath)
  {
    if (!cbica::directoryExists(dirName))
    {
      std::cerr << "Supplied directory name wasn't found: " << dirName << std::endl;
      exit(EXIT_FAILURE);
    }
    std::vector< std::string > allDirectories;
    std::string dirName_wrap = cbica::normPath(dirName);
    if (dirName_wrap[dirName_wrap.length() - 1] != '/')
    {
      dirName_wrap.append("/");
    }
#if defined(_WIN32)
    dirName_wrap.append("*.*");
    char* search_path = cbica::constCharToChar(dirName_wrap.c_str());
    WIN32_FIND_DATA fd;
    HANDLE hFind = ::FindFirstFile(search_path, &fd);
    if (hFind != INVALID_HANDLE_VALUE)
    {
      do
      {
        if ((fd.dwFileAttributes | FILE_ATTRIBUTE_DIRECTORY) == FILE_ATTRIBUTE_DIRECTORY && (fd.cFileName[0] != '.'))
        {
          if (returnFullPath)
          {
            allDirectories.push_back(dirName + "/" + std::string(fd.cFileName));
          }
          else
          {
            allDirectories.push_back(std::string(fd.cFileName));
          }
          if (recursiveSearch)
          {
            std::vector<std::string> tempVector = subdirectoriesInDirectory(dirName + "/" + std::string(fd.cFileName), true);
            allDirectories.insert(allDirectories.end(), tempVector.begin(), tempVector.end());
          }
        }
      } while (FindNextFile(hFind, &fd) != 0);
      ::FindClose(hFind);
    }
#else
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(dirName.c_str())) == NULL)
    {
      std::cerr << "Error(" << errno << ") occurred while opening directory '" << dirName << "'\n";
    }

    while ((dirp = readdir(dp)) != NULL)
    {
      if (recursiveSearch && (dirp->d_type == DT_DIR) && (dirp->d_name[0] != '.') && (dirp->d_name != std::string(".svn").c_str()))
      {
        std::vector<std::string> tempVector = subdirectoriesInDirectory(dirName + "/" + dirp->d_name, true);
        allDirectories.insert(allDirectories.end(), tempVector.begin(), tempVector.end());
      }

      if ((strcmp(dirp->d_name, ".") == 0) || (strcmp(dirp->d_name, "..") == 0))
        continue;

      if (dirp->d_type == DT_DIR)
      {
        allDirectories.push_back(dirName + "/" + dirp->d_name);
        if (returnFullPath)
        {
          allDirectories.push_back(dirName + "/" + dirp->d_name);
        }
        else
        {
          allDirectories.push_back(dirp->d_name);
        }
      }
    }
    closedir(dp);
#endif
    return allDirectories;
  }

  size_t numberOfRowsInFile(const std::string &csvFileName, const std::string &delim)
  {
    std::ifstream inFile(csvFileName.c_str());

    // new lines will be skipped
    inFile.unsetf(std::ios_base::skipws);

    // count the "\n"s with an algorithm specialized for counting
    return std::count(std::istream_iterator<char>(inFile), std::istream_iterator<char>(), *constCharToChar(delim));
  }

  size_t numberOfColsInFile(const std::string &csvFileName, const std::string & delim)
  {
    std::vector< std::string > rowVec;
    std::ifstream inFile(csvFileName.c_str());
    std::string line;
    std::getline(inFile, line, '\n');
    line.erase(std::remove(line.begin(), line.end(), '"'), line.end());

    // read a single row
    rowVec = stringSplit(line, delim);

    return rowVec.size();
  }

  std::vector< CSVDict > parseCSVFile(const std::string &csvFileName, const std::string &inputColumns, const std::string &inputLabels,
    bool checkFile, bool pathsRelativeToCSV, const std::string &rowsDelimiter, const std::string &colsDelimiter, const std::string &optionsDelimiter)
  {
    // if CSV file doesn't exist, exit with meaningful message
    if (!cbica::fileExists(csvFileName))
    {
      std::cerr << "Supplied file name, '" << csvFileName << "' wasn't found.\n";
      exit(EXIT_FAILURE);
    }

    //std::string dirName_Wrap = dirName;

    // store number of rows in the file - this is used to make the program parallelize-able
    const size_t numberOfRows = numberOfRowsInFile(csvFileName);

    // initialize return dictionary
    std::vector< CSVDict > return_CSVDict;
    return_CSVDict.resize(numberOfRows - 1);

    std::vector< std::vector < std::string > > allRows; // store the entire data of the CSV file as a vector of colums and rows (vector< rows <cols> >)

    std::ifstream inFile(csvFileName.c_str());

    for (size_t i = 0; i < numberOfRows; i++)
    {
      std::string line;
      std::getline(inFile, line, *constCharToChar(rowsDelimiter));
      line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
      allRows.push_back(stringSplit(line, colsDelimiter));
    }

    inFile.close(); // at this point, the entire data from the CSV file has been read and stored in allRows

    std::vector< std::string > inputColumnsVec = stringSplit(inputColumns, optionsDelimiter), inputLabelsVec; // columns to consider as images

    std::vector< size_t > inputColumnIndeces, // indeces in the CSV file which correspond to images
      inputLabelIndeces; // indeces in the CSV file which correspond to labels

    inputColumnIndeces.resize(inputColumnsVec.size());

    bool allLabelsConsidered = inputLabels.empty();
    if (!allLabelsConsidered)
    {
      inputLabelsVec = stringSplit(inputLabels, optionsDelimiter); // columns to consider as labels
      // if label columns are defined, use that number for the vector, otherwise, it considers all columns appearing AFTER the last image column as labels
      inputLabelIndeces.resize(inputLabelsVec.size());
    }

    // populate information about which indeces to store for data from first row (rowCounter = 0)
    for (size_t j = 0; j < inputColumnsVec.size(); j++)
    {
      inputColumnIndeces[j] = std::find(allRows[0].begin(), allRows[0].end(), inputColumnsVec[j]) - allRows[0].begin();
    }

    // populate information about which indeces to store for labels from first row (rowCounter = 0)
    if (!allLabelsConsidered)
    {
      // if specific labels are defined, use those indeces only
      for (size_t j = 0; j < inputLabelsVec.size(); j++)
      {
        inputLabelIndeces[j] = std::find(allRows[0].begin(), allRows[0].end(), inputLabelsVec[j]) - allRows[0].begin();
#ifdef __GNUC__
        inputLabelIndeces[j] -= inputLabelIndeces[j]; // this is done because gcc, for some weird reason, sets the absolute value for the indeces
#endif
      }
    }
    else
    {
      // otherwise, assume ALL other columns appearing AFTER the last image as labels
      for (size_t i = inputLabelIndeces.size() - 1; i < allRows[0].size(); i++)
      {
        inputLabelIndeces.push_back(i);
      }
    }

#ifdef _OPENMP
    // organize the data
    int threads = omp_get_max_threads(); // obtain maximum number of threads available on machine
    // if the total number of rows in CSV file are less than the available number of threads on machine (happens for testing),
    // use only the number of rows where meaningful data is present - this avoids extra thread overhead
    //threads > static_cast<int>(numberOfRows) ? threads = static_cast<int>(numberOfRows - 1) : threads = threads;
    //#pragma omp parallel for num_threads(threads)
#endif
    for (int rowCounter = 1; rowCounter < static_cast<int>(allRows.size()); rowCounter++)
    {
      return_CSVDict[rowCounter - 1].inputImages.resize(inputColumnIndeces.size()); // pre-initialize size to ensure thread-safety
      return_CSVDict[rowCounter - 1].inputLabels.resize(inputLabelIndeces.size()); // pre-initialize size to ensure thread-safety

      if (inputLabelIndeces.size() == 0) // contingency in the case where no labels are detected - just push back 1s instead
      {
        return_CSVDict[rowCounter - 1].inputLabels.resize(1);
      }
      for (size_t i = 0; i < inputColumnIndeces.size(); i++)
      {
        std::string fileToAdd = allRows[rowCounter][inputColumnIndeces[i]]; // case where the file names in the CSV are complete paths

        if (!fileExists(fileToAdd)) // if absolute paths aren't found, check for other conditions
        {
          if (pathsRelativeToCSV) // image paths are relative to location of CSV file
          {
            fileToAdd = cbica::getFilenamePath(csvFileName) + allRows[rowCounter][inputColumnIndeces[i]];
          }
          //else if (fileExists(dirName + fileToAdd))
          //{
          //  fileToAdd = dirName + fileToAdd;
          //}
          else // image paths are relative to CWD
          {
            fileToAdd = cbica::getCWD() + allRows[rowCounter][inputColumnIndeces[i]];
          }
        }

        return_CSVDict[rowCounter - 1].inputImages[i] = fileToAdd;

        if (checkFile) // this case should only be used for testing purposes
        {
          if (!fileExists(fileToAdd))
          {
            std::cerr << "File name in list does not exist. Location: row = " << rowCounter << ", col = " << inputColumnIndeces[i] << "\n";
            exit(EXIT_FAILURE);
          }
        }

        if (!inputLabelIndeces.empty())
        {
          for (size_t j = 0; j < inputLabelIndeces.size(); j++)
          {
            return_CSVDict[rowCounter - 1].inputLabels[j] = std::atof(allRows[rowCounter][inputLabelIndeces[j]].c_str());
          }
        }
        else
        {
          return_CSVDict[rowCounter - 1].inputLabels[0] = 1;
        }
      }

    }

    return return_CSVDict;
  }

  std::vector< std::vector< std::string > > readCSVDataFile(const std::string &csvFileName)
  {
    const size_t rows = numberOfRowsInFile(csvFileName);
    const size_t cols = numberOfColsInFile(csvFileName);

    std::vector< std::vector< std::string > > returnVector;
    std::ifstream data(csvFileName.c_str());
    std::string line, cell;

    returnVector.resize(rows);
    size_t i = 0, j = 0;
    while (std::getline(data, line))
    {
      j = 0;
      returnVector[i].resize(cols);
      std::stringstream lineStream(line);
      while (std::getline(lineStream, cell, ','))
      {
        returnVector[i][j] = cell;
        j++;
      }
      i++;
    }

    return returnVector;
  }

  std::map< std::string, size_t > ConfusionMatrix(const std::vector< float > &inputRealLabels, const std::vector< float > &inputPredictedLabels)
  {
    std::map< std::string, size_t > returnConfusionMatrix;

    if (inputRealLabels.size() != inputPredictedLabels.size())
    {
      std::cerr << "The sizes of the real and predicted labels do not match; exiting.\n";
      return returnConfusionMatrix;
    }

    size_t TP = 0, TN = 0, FP = 0, FN = 0, RP = 0, PP = 0;

    for (size_t i = 0; i < inputRealLabels.size(); i++)
    {
      if (inputRealLabels[i] == 1)
      {
        RP++;
      }
      if (inputPredictedLabels[i] == 1)
      {
        PP++;
      }

      // both real and predicted labels are equal means it is a "true" prediction
      if (inputRealLabels[i] == inputPredictedLabels[i])
      {
        if (inputRealLabels[i] == 1)
        {
          TP++;
        }
        else
        {
          TN++;
        }
      }
      else
      {
        if (inputRealLabels[i] == 1)
        {
          FN++;
        }
        else
        {
          FP++;
        }
      }
    }

    // construct the return structure
    returnConfusionMatrix["TP"] = TP;
    returnConfusionMatrix["FP"] = FP;
    returnConfusionMatrix["TN"] = TN;
    returnConfusionMatrix["FN"] = FN;
    returnConfusionMatrix["RP"] = RP;
    returnConfusionMatrix["PP"] = PP;

    return returnConfusionMatrix;
  }

  std::map< std::string, float > ROC_Values(const std::vector< float > &inputRealLabels, const std::vector< float > &inputPredictedLabels)
  {
    std::map< std::string, float > returnStatistics;

    if (inputRealLabels.size() != inputPredictedLabels.size())
    {
      std::cerr << "The sizes of the real and predicted labels do not match; exiting.\n";
      return returnStatistics;
    }

    auto confusionMatrix = ConfusionMatrix(inputRealLabels, inputPredictedLabels);
    returnStatistics["TP"] = static_cast<float>(confusionMatrix["TP"]);
    returnStatistics["FP"] = static_cast<float>(confusionMatrix["FP"]);
    returnStatistics["TN"] = static_cast<float>(confusionMatrix["TN"]);
    returnStatistics["FN"] = static_cast<float>(confusionMatrix["FN"]);
    returnStatistics["RP"] = static_cast<float>(confusionMatrix["RP"]);
    returnStatistics["PP"] = static_cast<float>(confusionMatrix["PP"]);

    // https://en.wikipedia.org/wiki/Accuracy_and_precision
    returnStatistics["Accuracy"] = (returnStatistics["TP"] + returnStatistics["TN"]) / (2 * inputRealLabels.size());

    // https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
    returnStatistics["PPV"] = returnStatistics["TP"] / returnStatistics["PP"];
    returnStatistics["Precision"] = returnStatistics["PPV"];

    // https://en.wikipedia.org/wiki/False_discovery_rate
    returnStatistics["FDR"] = returnStatistics["TP"] / returnStatistics["PP"];

    // https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values#false_omission_rate
    returnStatistics["FOR"] = returnStatistics["FN"] / (inputPredictedLabels.size() - returnStatistics["PP"]);

    // https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values
    returnStatistics["NPV"] = returnStatistics["TN"] / (inputPredictedLabels.size() - returnStatistics["PP"]);

    // https://en.wikipedia.org/wiki/Prevalence
    returnStatistics["Prevalence"] = returnStatistics["RP"] / (2 * inputRealLabels.size());

    // https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    returnStatistics["TPR"] = returnStatistics["TP"] / returnStatistics["RP"];
    returnStatistics["Sensitivity"] = returnStatistics["TPR"];
    returnStatistics["Recall"] = returnStatistics["TPR"];
    returnStatistics["POD"] = returnStatistics["TPR"];

    // https://en.wikipedia.org/wiki/False_positive_rate
    returnStatistics["FPR"] = returnStatistics["FP"] / (inputPredictedLabels.size() - returnStatistics["RP"]);
    returnStatistics["Fall-Out"] = returnStatistics["FPR"];

    // https://en.wikipedia.org/wiki/False_positives_and_false_negatives#False_positive_and_false_negative_rates
    returnStatistics["FNR"] = returnStatistics["FN"] / returnStatistics["RP"];
    returnStatistics["MR"] = returnStatistics["FNR"];

    // https://en.wikipedia.org/wiki/Sensitivity_and_specificity
    returnStatistics["TNR"] = returnStatistics["TN"] / returnStatistics["RP"];
    returnStatistics["Specificity"] = returnStatistics["TNR"];

    // https://en.wikipedia.org/wiki/Likelihood_ratios_in_diagnostic_testing#positive_likelihood_ratio
    returnStatistics["LR+"] = returnStatistics["TPR"] / returnStatistics["FPR"];

    // https://en.wikipedia.org/wiki/Likelihood_ratios_in_diagnostic_testing#negative_likelihood_ratio
    returnStatistics["LR-"] = returnStatistics["FNR"] / returnStatistics["TNR"];

    // https://en.wikipedia.org/wiki/Diagnostic_odds_ratio
    returnStatistics["DOR"] = returnStatistics["LR+"] / returnStatistics["LR-"];

    // https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
    returnStatistics["Dice"] = 2 * returnStatistics["TP"] / (2 * returnStatistics["TP"] + returnStatistics["FP"] + returnStatistics["FN"]);

    // https://en.wikipedia.org/wiki/Jaccard_index
    returnStatistics["JR"] = 2 * returnStatistics["TP"] / (returnStatistics["TP"] + returnStatistics["FP"] + returnStatistics["FN"]);

    return returnStatistics;
  }

  //inline std::string iterateOverStringAndSeparators(const std::string &inputString, size_t &count, int enum_separator = 10)
  //{
  //  std::string returnString = "";
  //  if (enum_separator == 10) // to get description
  //  {
  //    returnString = inputString.substr(count + 1); // get all characters after the separator was detected
  //  }
  //  else // for everything other than description
  //  {
  //    char testChar = inputString[count], separatorChar = *cbica::constCharToChar(getSeparator(enum_separator));
  //    size_t position, // position to start getting the substring
  //      separatorChecker = 2; // the configuration file needs the difference between two types of strings to be a single space (apart from the separator string)
  //    if (testChar == separatorChar)
  //    {
  //      count++;
  //      position = count;
  //      //stringStream.clear();
  //      //stringStream << inputString[count];
  //      //std::string testStr = stringStream.str(), testSep = getSeparator(enum_separator);
  //      //while (stringStream.str() != getSeparator(enum_separator))
  //      //{
  //      //  stringStream << inputString[count];
  //      //  returnString += stringStream.str();
  //      //  stringStream.clear();
  //      //  count++;
  //      //}
  //      testChar = inputString[count];
  //      while (testChar != separatorChar)
  //      {
  //        count++;
  //        testChar = inputString[count];
  //      }
  //
  //      returnString = inputString.substr(position, count - position);
  //    }
  //    else // a small check as a contingency plan
  //    {
  //      while (separatorChecker > 0)
  //      {
  //        separatorChecker--;
  //        count++;
  //        testChar = inputString[count];
  //      }
  //    }
  //  }

  //  return returnString;
  //}

  std::string getCurrentLocalDate()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    localtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%d:%02d:%02d",
      timeinfo.tm_year + 1900, timeinfo.tm_mon + 1, timeinfo.tm_mday);
#else
    tm *time_struct = localtime(&timer);
    sprintf(buffer, "%d:%02d:%02d",
      time_struct->tm_year + 1900, time_struct->tm_mon + 1, time_struct->tm_mday);
#endif

    return buffer;
  }

  std::string getCurrentLocalTime()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    localtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%02d:%02d:%02d",
      timeinfo.tm_hour, timeinfo.tm_min, timeinfo.tm_sec);
#else
    tm *time_struct = localtime(&timer);
    sprintf(buffer, "%02d:%02d:%02d",
      time_struct->tm_hour, time_struct->tm_min, time_struct->tm_sec);
#endif

    return buffer;
  }

  std::string getCurrentLocalDateAndTime()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    localtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%d:%02d:%02d,%02d:%02d:%02d",
      timeinfo.tm_year + 1900, timeinfo.tm_mon + 1, timeinfo.tm_mday,
      timeinfo.tm_hour, timeinfo.tm_min, timeinfo.tm_sec);
#else
    tm *time_struct = localtime(&timer);
    sprintf(buffer, "%d:%02d:%02d,%02d:%02d:%02d",
      time_struct->tm_year + 1900, time_struct->tm_mon + 1, time_struct->tm_mday,
      time_struct->tm_hour, time_struct->tm_min, time_struct->tm_sec);
#endif

    return buffer;
  }

  std::string getCurrentGMTDate()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    gmtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%d:%02d:%02d",
      timeinfo.tm_year + 1900, timeinfo.tm_mon + 1, timeinfo.tm_mday);
#else
    tm *time_struct = gmtime(&timer);
    sprintf(buffer, "%d:%02d:%02d",
      time_struct->tm_year + 1900, time_struct->tm_mon + 1, time_struct->tm_mday);
#endif

    return buffer;
  }

  std::string getCurrentGMT()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    gmtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%02d:%02d:%02d",
      timeinfo.tm_hour, timeinfo.tm_min, timeinfo.tm_sec);
#else
    tm *time_struct = gmtime(&timer);
    sprintf(buffer, "%02d:%02d:%02d",
      time_struct->tm_hour, time_struct->tm_min, time_struct->tm_sec);
#endif

    return buffer;
  }

  std::string getCurrentGMTDateAndTime()
  {
    time_t timer;
    // obtain current time
    time(&timer);
    char buffer[200];

    // obtain current local date
#ifdef _WIN32
    struct tm timeinfo;
    gmtime_s(&timeinfo, &timer);
    sprintf_s(buffer, "%d:%02d:%02d,%02d:%02d:%02d",
      timeinfo.tm_year + 1900, timeinfo.tm_mon + 1, timeinfo.tm_mday,
      timeinfo.tm_hour, timeinfo.tm_min, timeinfo.tm_sec);
#else
    tm *time_struct = gmtime(&timer);
    sprintf(buffer, "%d:%02d:%02d,%02d:%02d:%02d",
      time_struct->tm_year + 1900, time_struct->tm_mon + 1, time_struct->tm_mday,
      time_struct->tm_hour, time_struct->tm_min, time_struct->tm_sec);
#endif

    return buffer;
  }

  std::string getCurrentYear()
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

  //====================================== String stuff ====================================//

  bool splitFileName(const std::string &dataFile, std::string &path,
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
        dataFile_wrap = cbica::replaceString(dataFile_wrap, compressionFormats[i], "");
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
      extension = std::string(ext);
#else
      char *basename_var, *path_name; 

      auto idx = dataFile_wrap.rfind('.');

      // fun fact, replaceString(string, "", "") enters an infinite loop.
      // can we not use std::replace? is that not an option?
      if (extension != "") 
        dataFile_wrap = replaceString(dataFile_wrap, extension, "");
      
      if (idx != std::string::npos)
      {
        extension = "." + dataFile_wrap.substr(idx + 1);
        if (extension.find("/") != std::string::npos)
          extension = "";
        else {
#ifndef __APPLE__
          dataFile_wrap = replaceString(dataFile_wrap, extension, "");
#endif
        }
      }
      // else // there is no extension for file

      path_name = dirname(cbica::constCharToChar(dataFile_wrap.c_str()));
      basename_var = basename(cbica::constCharToChar(dataFile_wrap.c_str()));
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
      path = cbica::replaceString(path, "\\", "/"); // normalize path for Windows

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

#if (_MSC_VER >= 1700)
      path_name[0] = NULL;
      basename_var[0] = NULL;
      drive_letter[0] = NULL;
#endif
      if (!path.empty())
      {
        if (path[path.length() - 1] != '/')
        {
          path += "/";
        }
      }

      return true;
    }
  }

  std::vector<std::string> stringSplit(const std::string &str, const std::string &delim)
  {
    std::vector<std::string> results;

    for (size_t i = 0; i < str.length(); i++)
    {
      std::string tempString = "";
      while ((str[i] != *delim.c_str()) && (i < str.length()))
      {
        tempString += str[i];
        i++;
      }
      results.push_back(tempString);
    }

    return results;
  }

  std::string replaceString(const std::string &entireString,
    const std::string &toReplace,
    const std::string &replaceWith)
  {
    std::string return_string = entireString;
    for (size_t pos = 0;; pos += replaceWith.length())
    {
      pos = return_string.find(toReplace, pos);
      // std::cout << "pos: " << pos << std::endl;
      if (pos == std::string::npos)
        break;

      return_string.erase(pos, toReplace.length());
      return_string.insert(pos, replaceWith);
    }
    return return_string;
    /*
    if( entireString.length() < toReplace.length() )
    std::cerr << "Length of string to search < length of string to replace. Please check.\n";

    return(return_string.replace(entireString.find(toReplace), toReplace.length(), replaceWith));
    */
  }

  char* constCharToChar(const std::string &input)
  {
    char *s = new char[input.size() + 1];
#ifdef _WIN32
    strcpy_s(s, input.size() + 1, input.c_str());
#else
    std::strcpy(s, input.c_str());
#endif
    return s;
  }

  char* constCharToChar(const char *input)
  {
    return cbica::constCharToChar(std::string(input));
  }

  void dos2unix(const std::string inputFile)
  {
#ifndef WIN32 // this function is not needed for Windows systems
    auto tempDir = createTmpDir();
    auto tempFile = tempDir + "tempFile.txt";

    std::ifstream in(inputFile.c_str());
    if (!in.is_open())
    {
      std::cerr << "Error: could not open '" << inputFile << "'\n";
      return;
    }
    std::ofstream out(tempFile.c_str());
    std::istreambuf_iterator<char> input(in), end;
    std::ostreambuf_iterator<char> output(out);

    std::remove_copy(input, end, output, '\r');
    out.close();

    std::remove(inputFile.c_str());
    cbica::copyFile(tempFile, inputFile);

    if (removeDirectoryRecursively(tempDir) != 0)
    {
      std::cerr << "There was an issue deleting the tempDir '" << tempDir << "'\n";
    }
#endif
    return;
  }

  size_t getTotalMemory()
  {
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__))
    /* Cygwin under Windows. ------------------------------------ */
    /* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
    MEMORYSTATUS status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatus(&status);
    return (size_t)status.dwTotalPhys;

#elif defined(_WIN32)
    /* Windows. ------------------------------------------------- */
    /* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (size_t)status.ullTotalPhys;

#elif (defined(__APPLE__) && defined(__MACH__))
    	/* UNIX variants. ------------------------------------------- */
      /* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM */

    #if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
      int mib[2];
      mib[0] = CTL_HW;
      #if defined(HW_MEMSIZE)
        mib[1] = HW_MEMSIZE;            /* OSX. --------------------- */
      #elif defined(HW_PHYSMEM64)
        mib[1] = HW_PHYSMEM64;          /* NetBSD, OpenBSD. --------- */
      #endif
        int64_t size = 0;               /* 64-bit */
        size_t len = sizeof( size );
        if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
          return (size_t)size;
        return 0L;			/* Failed? */
    
    #endif
  #endif
  }

  size_t getCurrentlyUsedMemory()
  {
#if WIN32
    
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    GlobalMemoryStatusEx(&memInfo);
    return memInfo.ullTotalPhys - memInfo.ullAvailPageFile;

#elif (defined(__APPLE__) && defined(__MACH__))
    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics64_data_t vm_stats;

    mach_port = mach_host_self();
    count = sizeof(vm_stats) / sizeof(natural_t);
    if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
        KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
                                        (host_info64_t)&vm_stats, &count))
    {
        long long used_memory = ((int64_t)vm_stats.active_count +
                                 (int64_t)vm_stats.inactive_count +
                                 (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
        
        return (size_t)used_memory;
    }

#else
    return 0;
#endif
  }

  size_t getCurrentlyUsedMemoryByCurrentProcess()
  {
#if WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return pmc.WorkingSetSize;

#elif (defined(__APPLE__) && defined(__MACH__))
    
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL) 
    {
      if (strncmp(line, "VmRSS:", 6) == 0) 
      {
        result = strlen(line);
        const char* p = line;
        while (*p <'0' || *p > '9') 
          p++;
        line[result - 3] = '\0';
        result = atoi(p);
        break;
      }
    }
    fclose(file);
    return static_cast< size_t >(result * 1000);

#else
    return 0;
#endif
  }

  //! Cross platform Sleep
  void sleep(size_t ms)
  {
    if (ms <= 0)
    {
      std::cerr << "Sleep time of zero or less is not defined; calling Sleep with default parameter.\n";
      sleep();
    }

#if WIN32
    Sleep(static_cast<unsigned int>(ms));
#else
    #if __APPLE__
      time_t temp = static_cast<time_t>(ms);
    #else
      __time_t temp = static_cast<__time_t>(ms);
    #endif

    struct timespec ts = { temp / 1000, (temp % 1000) * 1000 * 1000 };
    nanosleep(&ts, NULL);
#endif
  }
}
