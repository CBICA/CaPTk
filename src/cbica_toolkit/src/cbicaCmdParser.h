/**
\file  cbicaCmdParser.h

\brief Declaration of the CmdParser class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include <string>
#include <stdio.h>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <map>

#if defined(__GNUC__) || defined(__clang__)
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

enum Separator
{
  Param, DataType, DataRange
};

//! String separators corresponding to Separator
#if defined(__GNUC__)  && (__GNUC__ < 5)
static const char *SeparatorStrings[] = { ":", "%", "*" };
#else
static std::vector< std::string > SeparatorStrings = { ":", "%", "*" };
#endif

//! Get the Separator as a string from the enum Separator
static inline std::string getSeparator(int enumVal)
{
  return SeparatorStrings[enumVal];
}

namespace cbica
{
  //! copied from cbicaUtilities to ensure CmdParser stays header-only
  inline std::string stringReplace(const std::string &entireString,
    const std::string &toReplace,
    const std::string &replaceWith)
  {
    std::string return_string = entireString;
    for (size_t pos = 0;; pos += replaceWith.length())
    {
      pos = return_string.find(toReplace, pos);
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

  //====================================== Structs that need string stuff ====================================//
  /**
  \struct Parameter

  \brief Holds individual parameter information

  This is a helper struct for internal usage of different functions and classes (right now, the function ReadConfigFile()
  and the class CmdParser() use it). It is not meant to be used from a program directly.
  All variables are self-explanatory. Currently, a maxium of five lines of description are supported.
  */
  struct Parameter
  {
    enum Type
    {
      FILE, DIRECTORY, STRING, INTEGER, FLOAT, BOOLEAN, NONE
    };

    std::string laconic;
    std::string verbose;
    int dataType_enumCode;
    std::string dataType_string;
    std::string dataRange;
    std::string descriptionLine1;
    std::string descriptionLine2; //! defaults to blank
    std::string descriptionLine3; //! defaults to blank
    std::string descriptionLine4; //! defaults to blank
    std::string descriptionLine5; //! defaults to blank

    size_t length;

    //! Constructor with five lines of description and enum_code for dataType
    Parameter(const std::string &in_laconic, const std::string &in_verbose, const int &in_dataType, const std::string &in_dataRange,
      const std::string &in_descriptionLine1, const std::string &in_descriptionLine2 = "", const std::string &in_descriptionLine3 = "",
      const std::string &in_descriptionLine4 = "", const std::string &in_descriptionLine5 = "") :
      laconic(in_laconic), verbose(in_verbose), dataType_enumCode(in_dataType), dataType_string(""), dataRange(in_dataRange),
      descriptionLine1(in_descriptionLine1), descriptionLine2(in_descriptionLine2),
      descriptionLine3(in_descriptionLine3), descriptionLine4(in_descriptionLine4), descriptionLine5(in_descriptionLine5)
    {
      laconic = cbica::stringReplace(laconic, "-", "");
      laconic = cbica::stringReplace(laconic, "--", "");
      verbose = cbica::stringReplace(verbose, "-", "");
      verbose = cbica::stringReplace(verbose, "--", "");
      length = laconic.length() + verbose.length();

      // populate dataType_string WRT dataType_enumCode
      switch (in_dataType)
      {
      case FILE:
        dataType_string = "FILE";
        break;
      case DIRECTORY:
        dataType_string = "DIRECTORY";
        break;
      case STRING:
        dataType_string = "STRING";
        break;
      case INTEGER:
        dataType_string = "INTEGER";
        break;
      case FLOAT:
        dataType_string = "FLOAT";
        break;
      case BOOLEAN:
        dataType_string = "BOOL";
        break;
      case NONE:
        dataType_string = "NONE";
        break;
      default:
        dataType_string = "UNKNOWN";
        break;
      }
    }

    //! Constructor with five lines of description and string for dataType
    Parameter(const std::string &in_laconic, const std::string &in_verbose, const std::string &in_dataType, const std::string &in_dataRange,
      const std::string &in_descriptionLine1, const std::string &in_descriptionLine2 = "", const std::string &in_descriptionLine3 = "",
      const std::string &in_descriptionLine4 = "", const std::string &in_descriptionLine5 = "") :
      laconic(in_laconic), verbose(in_verbose), dataType_enumCode(0), dataType_string(in_dataType), dataRange(in_dataRange),
      descriptionLine1(in_descriptionLine1), descriptionLine2(in_descriptionLine2),
      descriptionLine3(in_descriptionLine3), descriptionLine4(in_descriptionLine4), descriptionLine5(in_descriptionLine5)
    {
      laconic = cbica::stringReplace(laconic, "-", "");
      laconic = cbica::stringReplace(laconic, "--", "");
      verbose = cbica::stringReplace(verbose, "-", "");
      verbose = cbica::stringReplace(verbose, "--", "");
      length = laconic.length() + verbose.length();

      // populate dataType_enumCode WRT dataType_string
      if (dataType_string == "FILE")
      {
        dataType_enumCode = FILE;
      }
      else if (dataType_string == "DIRECTORY")
      {
        dataType_enumCode = DIRECTORY;
      }
      else if (dataType_string == "STRING")
      {
        dataType_enumCode = STRING;
      }
      else if (dataType_string == "INTEGER")
      {
        dataType_enumCode = INTEGER;
      }
      else if (dataType_string == "FLOAT")
      {
        dataType_enumCode = FLOAT;
      }
      else if ((dataType_string == "BOOL") || (dataType_string == "BOOLEAN"))
      {
        dataType_enumCode = BOOLEAN;
      }
      else if (dataType_string == "NONE")
      {
        dataType_enumCode = NONE;
      }
      else
      {
        dataType_enumCode = -1;
      }
    }

  };

  /**
  \class CmdParser

  \brief Simple command line parsing

  This is a pure c++ implementation. Executable name and project version are picked up automatically
  from the main CMakeLists file. Only the executable name can be modified in this class.

  An example of usage is shown below:

  \verbatim
  cbica::CmdParser parser = cbica::CmdParser(argc, argv); // OR,
  //cbica::CmdParser parser = cbica::CmdParser(argc, argv, "exe_name"); // if a different exe_name is desired

  /// The parameters "u"/"usage", "h"/"help" and "v"/"version" are automatically added ///

  // add parameters to the variable
  parser.addOptionalParameter("m","marvel", cbica::Parameter::INTEGER, "1 to 10", "I like The Avengers");
  parser.addOptionalParameter("d", "dc", cbica::Parameter::FLOAT, "1.00 to 10.00", "I prefer the Justice League");
  parser.addRequiredParameter("p", "people", cbica::Parameter::STRING, "max length = 1024", "People are always required");

  // note that there should be no spaces for either of the parameter types; they will be removed.

  /// checks for required parameters are done internally.

  std::string peopleString;
  parser.getParameterValue("p", peopleString);

  int marvelValue = 5; // set default value
  parser.getParameterValue("m", marvelValue);

  float dcValue = 5.15; // set default value
  parser.getParameterValue("d", dcValue);

  doSomethingWithTheParameters( peopleString, marvelValue, dcValue );
  \endverbatim
  */
  class CmdParser
  {
  public:
    /**
    \brief The Constructor

    \param argc The "argc" from executable
    \param argv The "argv" from executable
    \param exe_name Name of the executable, defaults to picking up from cbica::getExecutableName()
    */
    explicit CmdParser(const int argc, char **argv, const std::string &exe_name = "");

    /**
    \brief The Constructor

    \param argc The "argc" from executable
    \param argv The "argv" from executable
    \param exe_name Name of the executable, defaults to picking up from cbica::getExecutableName()
    */
    explicit CmdParser(const int argc, const char **argv, const std::string &exe_name = "");

    /**
    \brief The Destructor
    */
    virtual ~CmdParser();

    /**
    \brief Set a custom executable name
    */
    void setExeName(const std::string exeName){ m_exeName = exeName; };

    /**
    \get the executable name
    */
    std::string getExeName() { return m_exeName; };

    /**
    \brief Adding parameters: defaults to optional parameters

    As a standard, neither the laconic nor verbose parameters should have any '-' in the constructor.

    \param laconic The laconic variant
    \param verbose The verbose variant
    \param expectedDataType The data type expected for this parameter
    \param dataRange The range of data expected for this parameter
    \param description_line1 The first line of description for parameter
    \param description_line2 The second line of description for parameter, defaults to a blank string
    \param description_line3 The third line of description for parameter, defaults to a blank string
    \param description_line4 The fourth line of description for parameter, defaults to a blank string
    \param description_line5 The fifth line of description for parameter, defaults to a blank string
    */
    void addParameter(const std::string &laconic, const std::string &verbose, const int &expectedDataType, const std::string &dataRange,
      const std::string &description_line1, const std::string &description_line2 = "", const std::string &description_line3 = "",
      const std::string &description_line4 = "", const std::string &description_line5 = "");

    /**
    \brief Adding Optional parameters

    As a standard, neither the laconic nor verbose parameters should have any '-' in the constructor.

    \param laconic The laconic variant
    \param verbose The verbose variant
    \param expectedDataType The data type expected for this parameter
    \param dataRange The range of data expected for this parameter
    \param description_line1 The first line of description for parameter
    \param description_line2 The second line of description for parameter, defaults to a blank string
    \param description_line3 The third line of description for parameter, defaults to a blank string
    \param description_line4 The fourth line of description for parameter, defaults to a blank string
    \param description_line5 The fifth line of description for parameter, defaults to a blank string
    */
    void addOptionalParameter(const std::string &laconic, const std::string &verbose, const int &expectedDataType, const std::string &dataRange,
      const std::string &description_line1, const std::string &description_line2 = "", const std::string &description_line3 = "",
      const std::string &description_line4 = "", const std::string &description_line5 = "");

    /**
    \brief Adding Required parameters

    As a standard, neither the laconic nor verbose parameters should have any '-' in the constructor.

    \param laconic The laconic variant
    \param verbose The verbose variant
    \param expectedDataType The data type expected for this parameter
    \param dataRange The range of data expected for this parameter
    \param description_line1 The first line of description for parameter
    \param description_line2 The second line of description for parameter, defaults to a blank string
    \param description_line3 The third line of description for parameter, defaults to a blank string
    \param description_line4 The fourth line of description for parameter, defaults to a blank string
    \param description_line5 The fifth line of description for parameter, defaults to a blank string
    */
    void addRequiredParameter(const std::string &laconic, const std::string &verbose, const int &expectedDataType, const std::string &dataRange,
      const std::string &description_line1, const std::string &description_line2 = "", const std::string &description_line3 = "",
      const std::string &description_line4 = "", const std::string &description_line5 = "");

    /**
    \brief Display the usage
    */
    void echoUsage();

    /**
    \brief Display verbose usage
    */
    void echoHelp();

    /**
    \brief Display the version details
    */
    void echoVersion();

    /**
    \brief Check parameters WITHOUT hyphens

    Checks for both laconic and verbose variants of the specified parameter.

    \param execParamToCheck Which parameter to check
    \param position Position of parameter in argv else -1
    \return True if parameter found else False
    */
    bool compareParameter(const std::string &execParamToCheck, int &position, bool automaticEcho = true);

    /**
    \brief Check parameters WITHOUT hyphens

    Checks for both laconic and verbose variants of the specified parameter. Can be used to see if the parameter is present or not.

    \param execParamToCheck Which parameter to check
    \return True if parameter found else False
    */
    bool compareParameter(const std::string &execParamToCheck, bool automaticEcho = true);

    /**
    \brief Check if supplied parameter is present in the argument list

    Checks for both laconic and verbose variants of the specified parameter. Uses compareParameter() internally.

    \param execParamToCheck Which parameter to check
    \return True if parameter found else False
    */
    bool isPresent(const std::string &execParamToCheck, bool automaticEcho = true);

    /**
    \brief Get the description analogous with the parameter

    Can search using both laconic and verbose parameters.

    \param parameter Parameter whose description is requested
    \param "NewLine Return with "\n" between description lines if true, defaults to space between lines
    \return Description of parameter
    */
    std::string getDescription(const std::string &execParamToCheck, bool NewLine);

    /**
    \brief Get the data type analogous with the parameter

    Can search using both laconic and verbose parameters.

    \param parameter Parameter whose description is requested
    \return Description of parameter as string
    */
    std::string getDataTypeAsString(const std::string &execParamToCheck);

    /**
    \brief Get the data type analogous with the parameter

    Can search using both laconic and verbose parameters.

    \param parameter Parameter whose description is requested
    \return Description of parameter as Enum Code (Parameter::Type)
    */
    int getDataTypeAsEnumCode(const std::string &execParamToCheck);

    /**
    \brief Write the parser configuration as a JSON file
    */
    //void writeJSONFile(const std::string & dirName);

    /**
    \brief Write the configuration file for the executable for use in the common GUI framework

    The generated config file is always named 'EXE_NAME.txt'.

    \param dirName The full path of the directory to save the file; defaults to directory specified in cbica::makeTempDir()
    */
    void writeConfigFile(const std::string &dirName = "");

    /**
    \brief Reads a pre-written configuration file using CmdParser::WriteConfigFile()

    \param inputConfigFile Full path to the configuration file which needs to be read
    \return Vector of the Parameter structure where laconic paramter is always empty for all variables
    */
    static std::vector< Parameter > readConfigFile(const std::string &inputConfigFile, bool getDescription = true);

    /**
    \brief Gives a brief example of how to use the executable

    This should not contain any references to the executable name (it is automatically picked up).
    It should start directly with the parameters to be put in.

    \param usageOfExe A string which would correspond to the command line usage AFTER the executable has been called
    */
    DEPRECATED void exampleUsage(const std::string &usageOfExe);

    /**
    \brief Gives a brief example of how to use the executable

    This should not contain any references to the executable name (it is automatically picked up).
    It should start directly with the parameters to be put in.

    \param commandExcludingExeName A string which would correspond to the command line usage AFTER the executable has been called
    \param descriptionOfCommand A string which would correspond to what the command is expected to do
    */
    void addExampleUsage(const std::string &commandExcludingExeName, const std::string &descriptionOfCommand);

    /**
    \brief Adds descrition for the application for which the class is being initialized
    */
    void addApplicationDescription(const std::string &description);

	  /**
	  \brief Writes out a CWL specification file from cmd parser

	  This should be invoked everytime the cmd parser is called for an application.

	  \param dirName Full directory path to where the CWL spec will be produced
	  \param workflowName For more advanced CWL workflows
    \param overwriteFile If the current file should be overwritten or not
	  */
	  void writeCWLFile(const std::string & dirName, bool overwriteFile);

	  /*
	  \brief Gets the command line string with all the parameters from input yml file

	  This needs a input yml file to parse parameter values to CWL specs

	  \param dirName Full directory path to the input yml file
	
	  */
	  std::string GetCommandFromCWL(const std::string &inpDir, const std::string & cwlDir);

	  /*
	  \brief Reads a CWL spec file and creates a the command line tool for the application

	  This reads a CWL file and populates the cmd parser

	  \param path_to_config_file Path to the CWL spec file
	
	  \param getDescription This is a flag incase we want to run the CWL spec file with default values
	
	  */
	  void readCWLFile(const std::string & path_to_config_file, bool getDescription);

	  /*
	  \brief Creates a YAML node for the given root node.

	  This creates a new node

	  \param nodeString Name of the node

	  \param rootNode Parent node of the new node to be created

	  */
	  void createNode(const std::string & nodeString);

	  /*
	  \brief Deletes a YAML node inside the root node.

	  This creates a new node

	  \param nodeString Name of the node

	  \param parent Parent node of the new node to be created

	  */
	  void deleteNode(const std::string & nodeString);

	  /*
	  \brief Adds an input parameter inside input node.

	  This creates a new input parameter

	  \param param Name of the parameter to be added

	  \param rootNode Root Node of the CWL spec file

	  */
	  void addInputs(const std::string & param);

	  /*
	  \brief Adds an output field.

	  This creates a new output field

	  \param param Name of the parameter to be added

	  \param rootNode Root Node of the CWL spec file

	  */
	  void addOutputs(const std::string & param);

	  /*
	  \brief Checks if default value is specified for a parameter.

	  This returns the default value of parameter if exists

	  \param param Name of the parameter to be checked

	  \param input INPUT node of CWL spec

	  */
	  std::string checkDefault(const std::string & param);

	  /*
	  Invokes cwl-runner tool from the command line tool with given command

	  This invokes the cwl runner according to the flag for default values

	  \param cwl_spec_path Path to the CWL spec file
	  \param cwl_input_path Path to the input yml file. Empty string if doesn't exist
	  \param getDefaultFlag execute default or not default boolean

	  */
	  void cwlrunner(const std::string & cwl_spec_path, const std::string & cwl_input_path, bool getDefaultFlag);
	
	  /**
	  \brief Get the laconic value from verbose

	  Searches using the verbose value.

	  \param execParamToCheck The verbose variant of the parameter
	
	  */
	  std::string getLaconic(const std::string &execParamToCheck);

    /**
    \brief Get the value of the parameter

    Can search using both laconic and verbose parameters.

    \param execParamToCheck The laconic or verbose variant of the parameter
    \param parameterValue The return value of the parameter as bool
    */
    void getParameterValue(const std::string &execParamToCheck, bool &parameterValue);

    /**
    \brief Get the value of the parameter

    Can search using both laconic and verbose parameters.

    \param execParamToCheck The laconic or verbose variant of the parameter
    \param parameterValue The return value of the parameter as int
    */
    void getParameterValue(const std::string &execParamToCheck, int &parameterValue);

    /**
    \brief Get the value of the parameter

    Can search using both laconic and verbose parameters.

    \param execParamToCheck The laconic or verbose variant of the parameter
    \param parameterValue The return value of the parameter as size_t
    */
    void getParameterValue(const std::string &execParamToCheck, size_t &parameterValue);

    /**
    \brief Get the value of the parameter

    Can search using both laconic and verbose parameters.

    \param execParamToCheck The laconic or verbose variant of the parameter
    \param parameterValue The return value of the parameter as float
    */
    void getParameterValue(const std::string &execParamToCheck, float &parameterValue);

    /**
    \brief Get the value of the parameter

    Can search using both laconic and verbose parameters.

    \param execParamToCheck The laconic or verbose variant of the parameter
    \param parameterValue The return value of the parameter as std::string (valid for Parameter::Type::FILE, Parameter::Type::DIRECTORY, Parameter::Type::STRING)
    */
    void getParameterValue(const std::string &execParamToCheck, std::string &parameterValue);
    
    //! This function ensures that argc < 2 isn't checked
    void ignoreArgc1()
    {
      argc1ignore = true;
    }
    
  private:

    //! Executable name
    std::string m_exeName;
    //! Executable path
    std::string m_exePath;
    //! Version
    std::string m_version;
    //! Example of how to use the executable in question
    std::string m_exampleOfUsage;
    //! Description of the application - intended to give a short overview of what the user should expect
    std::string m_description;
    //! CMD variable, used to ensure that 'const' based variables are taken into consideration
    int m_argc;
    //! CMD variable, used to ensure that 'const' based variables are taken into consideration
    std::vector< std::string > m_argv;
    //! Collection of required and optional parameters
    std::vector< Parameter > m_requiredParameters, m_optionalParameters;
    //! Max length of parameters for echoUsage()
    size_t m_maxLength;
    //! Flag to toggle check for maximum overall length
    bool checkMaxLen;
    //! Flag to check for requested help/usage
    bool helpRequested;
    //! Flag to check for requested help/usage
    bool firstRun;
    //! check argc stuff internally
    bool argc1ignore;
    //! Initialize the class
    inline void initializeClass(int &input_argc, std::vector< std::string > &input_argv, const std::string &input_exeName = "");
    //! Get max length
    inline void getMaxLength();
    //! Internal function to check for verbose parameter
    inline void verbose_check(std::string &input_string);
    //! Internal function to write vector of parameters
    inline void writeParameters(const std::vector< Parameter > &inputParameters, bool verbose);
    //! Logs the CWL that is passed through cwl-runner
    void logCWL(const std::string &inpFileName, const std::string &cwlFileName);

    //! application usage examples and their respective descriptions
    std::vector< std::pair< std::string, std::string > > m_exampleUsageAndDescription;


    size_t m_maxLaconicLength, //! maximum length of laconic parameters
      m_minVerboseLength,
      m_maxVerboseLength; //! maximum length of verbose parameters

  };
}
