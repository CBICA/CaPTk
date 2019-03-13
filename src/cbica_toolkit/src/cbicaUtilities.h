/**
\file cbicaUtilities.h

\brief Some basic utility functions.

This needs c++11 flag enabled in gcc < 5.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include <string>
#include <typeinfo>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iterator>
#include <cmath>
#include <numeric>
#include <memory.h>
#include <map>
#include <random>
#include <iomanip>
#include <limits>

#if _WIN32
#include <process.h>
#endif

//#include <type_traits>

/**
\struct CSVDict

\brief A Dictionary structure to for CSV parsing

The CSV File in question is basically a list file consisting of file names (or paths) and labels
*/
struct CSVDict
{
  //! Contains input image file names
  std::vector< std::string > inputImages;

  //! Contains Labels that correspond to each file
  std::vector< double > inputLabels;

  //! Constructor
  CSVDict(const CSVDict &origin) :
    inputImages(origin.inputImages), inputLabels(origin.inputLabels)
  { };

  //! Constructor
  CSVDict()
  {
    inputImages.resize(0);
    inputLabels.resize(0);
  }

  //! Constructor
  CSVDict(const std::vector< std::string > &inputImagesVector, const std::vector< double > &inputLabelVector) :
    inputImages(inputImagesVector), inputLabels(inputLabelVector) {};
};

namespace cbica
{
  //====================================== String stuff ====================================//

  /**
  \brief Splits the input file name into its constituents

  \param dataFile The full file name which is the input
  \param baseName Overwritten with file name without extension
  \param extension Overwritten with extension without '.'
  \param path Overwritten with path to file

  \return True if successful
  */
  bool splitFileName(const std::string &dataFile, std::string &path,
    std::string &baseName, std::string &extension);

  /**
  \brief Splits the string

  \param str String to split
  \param delim Delimiter on the basis of which splitting is to be done
  \return results Output in the form of vector of strings
  */
  std::vector<std::string> stringSplit(const std::string &str, const std::string &delim);

  /**
  \brief Searches for smaller string in larger string and then replaces it with user-defined input

  \param entireString String to search
  \param toReplace String to replace
  \param replaceWith String to replace toReplace with

  \return std::string of result
  */
  std::string replaceString(const std::string &entireString,
    const std::string &toReplace,
    const std::string &replaceWith);

  /**
  \brief Convert const char* to char*

  \param input constant std::string

  \return character pointer
  */
  char* constCharToChar(const std::string &input);

  /**
  \brief Convert const char* to char*

  \param input constant character pointer

  \return character pointer
  */
  char* constCharToChar(const char *input);

  //======================================== OS stuff ======================================//

  /**
  \brief Check if file exists using istream

  \param fName Filename to check
  \return True if file exists
  */
  bool fileExists(const std::string &fName);

  /**
  \brief Check if directory exists

  \param dName String to check
  \return True if directory exists
  */
  bool directoryExists(const std::string &dName);

  /**
  \brief Return True if path is an existing regular file

  Reimplementation for python's "os.path.isfile": This follows symbolic links,
  so both islink() and isfile() can be true for the same path.

  [Wrap of cbica::fileExists()]

  \param path Filename of file to check

  \return True if path is an existing regular file
  */
  bool isFile(const std::string &path);

  /**
  \brief Return True if path is an existing directory

  Reimplementation of python's "os.path.isdir": This follows symbolic links,
  so both islink() and isdir() can be true for the same path.

  [Wrap of cbica::directoryExists()]

  \param path Directory name to check

  \return True of path is an existing directory
  */
  bool isDir(const std::string &path);

  /**
  \brief Return True if path exists and false for broken symbolic links

  Reimplementation of python's "os.path.exists": On some platforms, this function may return False
  if permission is not granted to execute os.stat() on the requested file, even if the
  path physically exists.

  \param path Path to check

  \return True if path is valid and false for broken symbolic links
  */
  bool exists(const std::string &path);

  std::vector<std::string> getCWLFilesInApplicationDir();

  /**
  \brief Create a temporary directory

  Creates a user-writable file using the following format: USER_HOME_DIR + EXE_NAME + tmp_ + processID;

  If this file is existing, for whatever reason, then a new is created which has the time stamp appended.

  \return Path of temporary directory
  */
  std::string createTmpDir();

  /**
  \brief Create a temporary directory

  Wrap for createTmpDir()

  \return True if success
  */
  std::string createTemporaryDirectory();

  /**
  \brief Create a temporary directory

  Wrap for createTmpDir()

  \return True if success
  */
  std::string makeTemporaryDirectory();

  /**
  \brief Create a temporary directory

  Wrap for createTmpDir()

  \return True if success
  */
  std::string makeTempDir();

  /**
  \brief Create a directory

  \param dir_name Name of directory to be created with full path

  \return True if success
  */
  bool createDir(const std::string &dir_name);

  /**
  \brief Create a directory

  \param dir_name Name of directory to be created with full path
  */
  bool makeDir(const std::string &dir_name);

  /**
  \brief Create a directory

  \param dir_name Name of directory to be created with full path

  \return True if success
  */
  bool createDirectory(const std::string &dir_name);

  /**
  \brief Create a directory

  \param dir_name Name of directory to be created with full path
  */
  bool makeDirectory(const std::string &dir_name);

  /**
  \brief Create a directory

  \param dir_name Name of directory to be created with full path

  \return True if success
  */
  bool createFolder(const std::string &dir_name);

  /**
  \brief Create a directory

  \param dir_name Name of directory to be created with full path
  */
  bool makeFolder(const std::string &dir_name);

  /**
  \brief Recursively delete a folder and contents [internal function]

  \param dirname Folder to delete

  \return true for success
  */
  int removeDirectoryRecursively(const std::string &dirname, bool bDeleteSubdirectories);

  /**
  \brief Delete a folder and contents

  \param path Folder to delete

  \return true for success
  */
  bool removeDir(const std::string &path);

  /**
  \brief Delete a folder and contents

  \param path Folder to delete

  \return true for success
  */
  bool deleteDir(const std::string &path);

  /**
  \brief Copy a folder and if recursion enabled, all its contents

  https://msdn.microsoft.com/en-us/library/hh874694.aspx?f=255&MSPPError=-2147217396

  \param inputFolder Folder to copy
  \param destination Where to copy to
  \param recursion Do recursion and copy, defaults to true

  \return true for success
  */
  bool copyDir(const std::string &inputFolder, const std::string &destination, bool recursion = true);

  /**
  \brief Copy a folder and if recursion enabled, all its contents

  https://msdn.microsoft.com/en-us/library/hh874694.aspx?f=255&MSPPError=-2147217396

  \param inputFolder Folder to copy
  \param destination Where to copy to
  \param recursion Do recursion and copy, defaults to true

  \return true for success
  */
  bool copyDirectory(const std::string &inputFolder, const std::string &destination, bool recursion = true);

  /**
  \brief Copy a folder and if recursion enabled, all its contents

  https://msdn.microsoft.com/en-us/library/hh874694.aspx?f=255&MSPPError=-2147217396

  \param inputFolder Folder to copy
  \param destination Where to copy to
  \param recursion Do recursion and copy, defaults to true

  \return true for success
  */
  bool copyFolder(const std::string &inputFolder, const std::string &destination, bool recursion = true);

  /**
  \brief Copy a folder and if recursion enabled, all its contents

  \param inputFile File to copy
  \param destination Where to copy to
  \return true for success
  */
  bool copyFile(const std::string &inputFile, const std::string &destination);

  /**
  \brief Get the size of the file in bytes

  \param inputFile The input file
  */
  size_t getFileSize(const std::string &inputFile);

  /*
  \brief Checks for the compatibility with the current project

  \param inputVersionFile The version file (in YAML) that contains the compatibility information 
  */
  bool IsCompatible(const std::string inputVersionFile);

  /**
  \brief Get the size of the folder

  \param rootFolder The input folder
  */
  size_t getFolderSize(const std::string &rootFolder);

  /**
  \brief Get the size of the folder

  \param rootFolder The input folder
  */
  size_t getDirSize(const std::string &rootFolder);

  /**
  \brief Get the size of the folder

  \param rootFolder The input folder
  */
  size_t getDirectorySize(const std::string &rootFolder);

  /**
  \brief Gets the extension of the supplied file name using splitFileName()

  Prefer to use "/" as file path delimiter.

  \param filename The input filename
  \param checkFile Checks existence of file using fileExists
  \return std::string which has the file extension
  */
  std::string getFilenameExtension(const std::string &filename, bool checkFile = true);

  /**
  \brief Gets the base of the supplied file name using splitFileName()

  Prefer to use "/" as file path delimiter.

  \param filename The input filename
  \param checkFile Checks existence of file using fileExists
  \return std::string which has the file extension
  */
  std::string getFilenameBase(const std::string &filename, bool checkFile = true);

  /**
  \brief Gets the path of the supplied file name using splitFileName()

  Prefer to use "/" as file path delimiter.

  \param filename The input filename
  \param checkFile Checks existence of file using fileExists
  \return std::string which has the file extension
  */
  std::string getFilenamePath(const std::string &filename, bool checkFile = true);

  /**
  \brief Get the name of the Executable which is calling the function

  \return exe name
  */
  std::string getExecutableName();

  /**
  \brief Get the path of the Executable which is calling the function

  returns the path after calling splitFileName on getFullPath

  \return exe name
  */
  std::string getExecutablePath();

  /**
  \brief Get the name of the Executable which is calling the function

  \return exe name
  */
  std::string getFullPath();

  /**
  \brief Get the name of the user who is calling the function

  \return user name
  */
  std::string getUserName();

  /**
  \brief Get the home directory of the user

  Windows - C:/Users/XYZ
  Linux - /home/XYZ

  \return user directory
  */
  std::string getUserHomeDirectory();

  /**
  \brief Get the current working directory

  [Wrap for cbica::getcwd()]

  \return Current working directory
  */
  std::string getCWD();

  /**
  \brief Normalize a pathname by collapsing redundant separators and up-level references

  Reimplementation of python's "os.path.normpath": Normalize a pathname by collapsing redundant separators
  and up-level references so that A//B, A/B/, A/./B and A/foo/../B all become A/B. This string
  manipulation may change the meaning of a path that contains symbolic links. On Windows, it
  converts forward slashes to backward slashes.

  \param path Path to normalize

  \return std::string Normalized path
  */
  std::string normPath(const std::string &path);

  /**
  \brief Normalize a pathname by collapsing redundant separators and up-level references

  [Wrap for cbica::normPath()]

  \param path Path to normalize

  \return std::string Normalized path
  */
  std::string normalizePath(const std::string &path);

  /**
  \brief Return a relative filepath to path

  Reimplementation for python's "os.path.relpath": Return a relative filepath to path either from the
  current directory or from an optional start directory. This is a path computation: the
  filesystem is not accessed to confirm the existence or nature of path or start.

  \param path Path to check
  \param base Base file name

  \return Relative filepath to path either from current directory or from an optional start dir
  */
  std::string relPath(const std::string &path, const std::string &base);

  /**
  \brief Return a relative filepath to path

  Wrap of cbica::relPath()

  \param path Path to check
  \param base Base file name

  \return Relative filepath to path either from current directory or from an optional start dir
  */
  std::string relativePath(const std::string &path, const std::string &base);

  /**
  \brief Return the canonical path of the specified filename

  Reimplementation of python's "os.path.realpath": Return the canonical path of the specified
  filename, eliminating any symbolic links encountered in the path (if they are
  supported by the operating system).

  \param path Filename to check

  \return canonical path of path
  */
  std::string realPath(const std::string &path);

  /**
  \brief Return True if path refers to a directory entry that is a symbolic link

  Reimplementation of python's "os.path.islink": Always False if symbolic links are not
  supported by the python runtime.

  \param path Path to check

  \return True if path is symbolic link
  */
  bool isLink(const std::string &path);

  /**
  \brief Check if path refers to a symbolic entry

  [Wrap of cbica::isLink()]

  \param path Path to check

  \return True if path is symbolic link
  */
  bool isSymbolicLink(const std::string &path);

  /**
  \brief Make a symbolic link from file to another

  \param input_fileName Input file name for which symbolic link needs to be created
  \param output_fileName Output file name which is the symbolic link for input_fileName

  \return True if symbolic link successfully created
  */
  bool makeSymbolicLink(const std::string &input_fileName, const std::string &ouput_fileName);

  /**
  \brief Sets the environment variable

  \param variable_name Name of the Variable
  \param variable_value Value of variable_name

  \return True if successful
  */
  bool setEnvironmentVariable(const std::string &variable_name, const std::string &variable_value);

  /**
  \brief Gets the value of the specified environment variable
  */
  std::string getEnvironmentVariableValue(const std::string &environmentVariable);

  /**
  \brief Delete the environment variable

  \param variable_name Name of the Variable

  \return True if successful
  */
  bool deleteEnvironmentVariable(const std::string &variable_name);

  /**
  \brief Find all files inside a directory

  \param dirName The directory to do the search in
  */
  std::vector< std::string > filesInDirectory(const std::string &dirName, bool returnFullPath = true);

  /**
  \brief Find all sub-directories inside a directory

  \param dirName The directory to do the search in
  \param recursiveSearch Whether to do a recursive search or on a single level
  */
  std::vector<std::string> subdirectoriesInDirectory(const std::string &dirName, bool recursiveSearch = false, bool returnFullPath = false);

  /**
  \brief Find number of rows in CSV file
  */
  size_t numberOfRowsInFile(const std::string &csvFileName, const std::string &delim = "\n");

  /**
  \brief Find number of cols in CSV file
  */
  size_t numberOfColsInFile(const std::string &csvFileName, const std::string &delim = ",");

  /**
  \brief Parse the supplied CSV File and obtain Row and Column information.

  Assumptions:
  1. Header information is in first row
  2. Full paths of images are given
  3. Image paths are either absolute (default behaviour) OR relative to the location of CSV file OR relative to CWD

  \param csvFileName The full path of the file to parse, all paths are absolute or relative to current working directory
  \param inputColumns The string of input columns which contain the data to be used for further processing
  \param inputLabels The string of input labels per subject based on which further processing is to be done; if this is empty, it is initialized as 1 for all subjects
  \param checkFile Check the validity of the file; defaults to true
  \param pathsRelativeToCSV The paths in the CSV file are relative to the location of CSV file otherwise, they are checked relative to the CWD
  \param rowsDelimiter The delimiters used to distinguish rows in the file
  \param colsDelimiter The delimiters used to distinguish cols in the file
  \param optionsDelimiter The delimiters used in inputColumns and inputLabel
  \return Vector of CSV Dictionary items: Collection of input images are respective labels
  */
  std::vector< CSVDict > parseCSVFile(const std::string &csvFileName, const std::string &inputColumns, const std::string &inputLabels, bool checkFile = true, bool pathsRelativeToCSV = false, const std::string &rowsDelimiter = "\n", const std::string &colsDelimiter = ",", const std::string &optionsDelimiter = ",");

  /**
  \brief Read a CSV file which has no header information

  To read CSV file with header information, check parseCSVFile() function. This should not be used for obtaining strings

  \param csvFileName The full path of the file to parse, all paths are absolute or relative to current working directory
  \param columnMajor If true, then return is a vector of all the columns; otherwise it is a vector of the rows
  */
  template< class TDataType = double >
  std::vector< std::vector< TDataType > > readCSVDataFile(const std::string &csvFileName, bool columnMajor = false)
  {
    std::vector< std::vector< TDataType > > returnVector;
    if (!cbica::isFile(csvFileName))
    {
      std::cerr << "Supplied file wasn't found.\n";
      return returnVector;
    }

    const size_t rows = numberOfRowsInFile(csvFileName);
    const size_t cols = numberOfColsInFile(csvFileName);

    std::ifstream data(csvFileName.c_str());
    std::string line, cell;

    if (columnMajor)
    {
      returnVector.resize(cols);
    }
    else
    {
      returnVector.resize(rows);
    }

    size_t i = 0, j = 0;
    while (std::getline(data, line))
    {
      j = 0;
      if (!columnMajor)
      {
        returnVector[i].resize(cols);
      }
      std::stringstream lineStream(line);
      while (std::getline(lineStream, cell, ','))
      {
        if (columnMajor && returnVector[j].empty())
        {
          returnVector[j].resize(rows);
        }

        auto temp = static_cast<TDataType>(std::atof(cell.c_str()));

        if (columnMajor)
        {
          returnVector[j][i] = temp;
        }
        else
        {
          returnVector[i][j] = temp;
        }
        j++;
      }
      i++;
    }

    return returnVector;
  }

  /**
  \brief Find the unique elements in a vector

  Implementation incorporated from the SO answer in https://stackoverflow.com/a/1041939/1228757

  \param inputVector The vector on which to search for unique elements 
  \return An std::vector with unique elements
  */
  template< class TDataType = std::string >
  std::vector< TDataType > GetUniqueElements(const std::vector< TDataType > &inputVector)
  {
    std::set< TDataType > s;
    std::vector< TDataType > returnVector;
    for (size_t i = 0; i < inputVector.size(); i++)
    {
      s.insert(inputVector[i]);
    }
    returnVector.assign(s.begin(), s.end());

    return returnVector;
  }

  /**
  \brief Read a CSV file which has no header information

  To read CSV file with header information, check parseCSVFile() function. This should be used for obtaining strings

  \param csvFileName The full path of the file to parse, all paths are absolute or relative to current working directory
  */
  std::vector< std::vector< std::string > > readCSVDataFile(const std::string &csvFileName);

  /**
  \brief Get current local time as string delineated as YYYY:MM:DD
  */
  std::string getCurrentLocalDate();

  /**
  \brief Get current local time as string delineated as HH:MM:SS
  */
  std::string getCurrentLocalTime();

  /**
  \brief Get current local time as string delineated as YYYY:MM:DD,HH:MM:SS
  */
  std::string getCurrentLocalDateAndTime();

  /**
  \brief Get current GMT as string delineated as YYYY:MM:DD
  */
  std::string getCurrentGMTDate();

  /**
  \brief Get current GMT as string delineated as HH:MM:SS
  */
  std::string getCurrentGMT();

  /**
  \brief Get current GMT as string delineated as YYYY:MM:DD,HH:MM:SS
  */
  std::string getCurrentGMTDateAndTime();

  /**
  \brief Get current Year as string delineated as YYYY
  */
  std::string getCurrentYear();

  /**
  \brief Get the current process ID

  Provides wraps to _getpid() in OS-specific ways
  */
  std::string getCurrentProcessID();

  /**
  \brief Cross platform sleep

  Defaults to "std::rand() % 1000 + 1"
  */
  void sleep(size_t ms = std::rand() % 1000 + 1);

  /**
  \brief Ensuring files written using Windows don't mess stuff up

  Base implementation from https://www.digitalpeer.com/blog/simple-text-processing-with-cpp-dos2unix-example
  */
  void dos2unix(const std::string inputFile);

  //==================================== Template stuff ==================================//

  /**
  \brief Searches for an element in a vector and returns true/false and position

  Templated function to take in any kind of vector and element.

  \param vector_to_search_in Vector to do the search in
  \param element_to_search_for Element to search for
  \param position Last position of found element in vector (-1) if not found

  \return True if found
  \return Position if found (-1) if not
  */
  template<typename TContainerType = std::string >
  std::pair<bool, int> findInVector(std::vector<TContainerType> &vector_to_search_in,
    TContainerType element_to_search_for)
  {
    int position = -1;
    //std::vector<int>::const_iterator
    typename std::vector<TContainerType>::const_iterator iterator =
      std::find(vector_to_search_in.begin(), vector_to_search_in.end(), element_to_search_for);
    if (iterator != vector_to_search_in.end())
    {
      position = iterator - vector_to_search_in.begin();
      return std::make_pair(true, position);
    }
    else
      return std::make_pair(false, position);
  }

  /**
  \brief Convert first character to integer, double, unsigned int, etc.

  \param input_string Input character to be converted
  \return Templated to the type of return required
  */
  template < typename TConvertType
#if (_MSC_VER >= 1800) || (__GNUC__ > 4)
    = int
#endif
>
/*typename*/ TConvertType convertCharacter(const std::string &input_string)
  {
    return static_cast<TConvertType>(input_string[0]);
  }

  /**
  \brief Convert entire string to integer, double, unsigned int, etc.

  \param input_string Input character to be converted
  \return Templated vector to the type of return required
  */
  template<typename TConvertType = int >
  std::vector</*typename*/ TConvertType> convertString(const std::string &input_string)
  {
    std::vector</*typename*/TConvertType>return_vector;
    for (int i = 0; i < input_string.length(); ++i)
      return_vector.push_back(static_cast<TConvertType>(input_string[i]));

    return return_vector;
  }

#if (_MSC_VER >= 1800) || (__cplusplus >= 201103L)
  /**
  \brief Base for compareEqual(...)
  */
  template <typename A = int, typename B = A>
  inline bool compareEqual(const A x, const B y)
  {
    return (x == y);
  }

  /**
  \brief Compare if equal for multiple inputs
  */
  template <typename A = int, typename B = A, typename... Others>
  bool compareEqual(const A x, const B y, Others const ... args)
  {
    return (x == y) && compareEqual(y, args...);
  }

  /**
  \brief Base for compareGreater(...)
  */
  template <typename A = int, typename B = A>
  inline bool compareGreater(const A x, const B y)
  {
    return (x > y);
  }

  /**
  \brief Compare if greater for multiple inputs
  */
  template <typename A = int, typename B = A, typename... Others>
  bool compareGreater(const A x, const B y, Others const ... args)
  {
    return (x > y) && compareGreater(y, args...);
  }

  /**
  \brief Base for compareLesser(...)
  */
  template <typename A = int, typename B = A>
  inline bool compareLesser(const A x, const B y)
  {
    return (x < y);
  }

  /**
  \brief Compare if lesser for multiple inputs
  */
  template <typename A = int, typename B = A, typename... Others>
  bool compareLesser(const A x, const B y, Others const ... args)
  {
    return (x < y) && compareLesser(y, args...);
  }

  /**
  \brief Convert to string with precision

  Useful when std::to_string truncates doubles
  */
  template< typename TDataType = double >
  std::string to_string_precision(const TDataType a_value, const int n = 10)
  {
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
  }

#endif

  //==================================== Statistical/Compute stuff ==================================//
  /**
  \brief Calculates the Confusion Matrix for a set of real and predicted labels

  Values returned: True Positive (TP), False Positive (FP), True Negative (TN), False Negative (FN), Real Positive (RP), Preditcted Positive (PP)

  \param inputRealLabels Vector structure containing real labels
  \param inputPredictedLabels Vector structure containing predicted labels

  \return std::map< string, size_t > A map of string and corresponding non-negative values
  */
  std::map< std::string, size_t > ConfusionMatrix(const std::vector< float > &inputRealLabels, const std::vector< float > &inputPredictedLabels);

  /**
  \brief Calculates the ROC Values (see https://en.wikipedia.org/wiki/Receiver_operating_characteristic of all estimates) for a set of real and predicted labels

  Values returned:
  True Positive (TP), False Positive (FP), True Negative (TN), False Negative (FN), Real Positive (RP), Preditcted Positive (PP) [all from ConfusionMatrix()],
  Accuracy, Positive Predictive Value (PPV) [aka Precision], False Discovery Rate (FDR), False Omission Rate (FOR),
  Negative Predictive Value (NPV), Prevalence, True Positive Rate (TPR) [aka Sensitivity, Recall, Probability of Detection (POD)],
  False Positive Rate (FPR) [aka Fall-out], False Negative Rate (FNR) [aka Miss Rate (MR)],  True Negative Rate (TNR) [aka Specificity],
  Positive Likelihood Ratio (LR+), Negative Likelihood Ratio (LR−), Diagnostic Odds Ratio (DOR), Dice Score (Dice), Jaccard ratio/index (JR)

  \param inputRealLabels Vector structure containing real labels
  \param inputPredictedLabels Vector structure containing predicted labels

  \return std::map< string, size_t > A map of string and corresponding float values
  */
  std::map< std::string, float > ROC_Values(const std::vector< float > &inputRealLabels, const std::vector< float > &inputPredictedLabels);

  /**
  \brief A good random number generator using c++11 that gives a random value within a range

  \param start Start value of range; defaults to 0
  \param end End value of range; defaults to 1
  \param sizeOfReturn Size of the returned vector
  */
  template< class TDataType = float >
  std::vector< TDataType > randn(const TDataType start, const TDataType end, size_t sizeOfReturn = 1)
  {
    std::vector< TDataType > returnVec;
    returnVec.resize(sizeOfReturn);
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::normal_distribution< TDataType > distr(start, end); // define the range

    for (size_t i = 0; i < returnVec.size(); i++)
    {
      returnVec[i] = distr(eng);
    }

    return returnVec;
  }

  /**
  \brief A good random number generator using c++11 that gives a random value within a range

  \param sizeOfReturn Size of the returned vector
  */
  template< class TDataType = float >
  std::vector< TDataType > randn(size_t sizeOfReturn = 1)
  {
    std::vector< TDataType > returnVec;
    returnVec.resize(sizeOfReturn);
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::normal_distribution< TDataType > distr(); // define the range

    for (size_t i = 0; i < returnVec.size(); i++)
    {
      returnVec[i] = distr(eng);
    }

    return returnVec;
  }

};

/**
\struct FileNameParts

\brief Holds the different parts of a file name (path, base and extension)

This is a helper struct for uses concerning uses where different parts of a file (path, base, ext) need to be used multiple times.
*/
struct FileNameParts
{
  std::string fullFileName, path, base, extension;

  //! Default Constructor
  FileNameParts();

  /**
  \brief Constructor with input file string

  \param inputFileName The complete input file name
  */
  FileNameParts(const std::string &inputFileName) : fullFileName(inputFileName)
  {
    cbica::replaceString(fullFileName, "\\", "/");
    if (cbica::fileExists(inputFileName))
      cbica::splitFileName(fullFileName, path, base, extension);
    else
    {
      std::cerr << "The input file '" << inputFileName << "' wasn't found on disk.\n";
    }
  }

  /**
  \brief Member function to set the fullFileName

  \param inputFileName The complete input file name
  */
  void SetFileName(const std::string &inputFileName)
  {
    fullFileName = inputFileName;
    cbica::replaceString(fullFileName, "\\", "/");
    if (cbica::fileExists(fullFileName))
      cbica::splitFileName(fullFileName, path, base, extension);
    else
    {
      std::cerr << "The input file '" << inputFileName << "' wasn't found on disk.\n";
    }
  }
};
