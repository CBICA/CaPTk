/**
\file  cbicaLogging.h

\brief Declaration of the Logging class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include <iostream>
#include <fstream>
#include <memory>
#include <time.h>
#include <stdio.h>

#include "cbicaUtilities.h"

/*
\namespace cbica
\brief Namespace for differentiating functions written for internal use
*/
namespace cbica
{
  /**
  \class Logging

  \brief The logging class.

  This automatically generates a machine-parseable log specified by the file name. The user also has the option 
  of submitting free text to be put along with the log. The generated log is in the format show below:

  <code><4 digit year>:<2 digit month>:<2 digit date>,<2 digit 24 hour>:<2 digit minute>:<2 digit second>;<exe name>;<user name>;<free text></code>

  Usage:
  \verbatim
  // writing to a file
  cbica::Logging logger( "file_name.txt", "randomness is highly underrated" ); // the file has already been written at this point and can be viewed
  // to write to console, either initialize the class as 'cbica::Logging logger;' or call EnableConsoleLogging() after the initialization.
  ...
  // more code
  ...
  logger.Write( "'I accept chaos, I'm not sure whether it accepts me' - Bob Dylan" ); // writes to file_name.txt
  \endverbatim

  The class defaults to console logging. Use Logging::EnableTextLogging() or UseNewFile() to switch; 
  by default, it writes to a text file 'EXE_NAME-log.txt' in directory specified by cbica::createTemporaryDirectory()
  */
	class Logging
	{
	public:
		/**
    \brief Actual Constructor

    \param file_name_with_path The file onto which the log file is to be written
    \param FreeText_input Free text which the user wants to be present in the log, defaults to an empty string
    */
		explicit Logging(const std::string file_name, const std::string FreeText_input);
		
		/**
    \brief Default constructor

    Just used to keep a track of the user name and executable run at a particular time.
    */
		explicit Logging();

    /**
    \brief Default constructor

    Just used to keep a track of the user name and executable run at a particular time.
    */
    Logging(const Logging &origin);

    /**
    \brief Change Logging file after initializing class

    \param newLogFile Path of new log file. If empty, it becomes 'cbica::createTmpDir() + cbica::getExecutableName() + "-log.txt"'
    */
    void UseNewFile(const std::string &newLogFile);

		//! The Destructor
		virtual ~Logging();

    /**
    \brief Function to call to write error messages to log file without any free text

    \param FreeText_input Free text which the user wants to be present in the log
    */
    void WriteError(const std::string FreeText_input);

    /**
		\brief Function to call to write to log file
        
    \param FreeText_input Free text which the user wants to be present in the log, defaults to an empty string
    */
		void Write(const std::string FreeText_input);
    
    /**
    \brief Switches from console to text file logging
    
    The output stamps are of the form:

    <4 digit year>:<2 digit month>:<2 digit date>,<2 digit 24 hour>:<2 digit minute>:<2 digit second>;<free text>
    */
    void EnableTextLogging(const std::string &newLogFile);
    
    /**
    \brief This is useful in scenarios where the application has been installed in a multi-user machine.

    This also disables writing of the date in the log since the assumption is that on a single user machine, the filename should be in the following format: 
    ${userHomeDir}/.${appName}/${currentDate} (see CaPTk logging as example). This also disables the writing of the executable name.
    */
    void EnableMultiUserLogging()
    {
      m_multiUserLog = true;
    }

    /**
    \brief Switches from text to console file logging

    This is helpful if the user wants to visualize the console output. If it is done for saving, the recommended way is to call EnableTextLogging().

    The output stamps are of the form:

    <4 digit year>:<2 digit month>:<2 digit date>,<2 digit 24 hour>:<2 digit minute>:<2 digit second>;<free text>
    */
    void EnableConsoleLogging();

    /**
    \brief This enables logging the date and time in GMT rather than in local (which is the default behavior)
    */
    void EnableGMTLogging();

    /**
    \brief Get the file name with full path where log has happened

    \return file_name_with_path 
    */
    std::string getLoggingFileName();

  protected:    
    /**
    \brief The function used to initialize the class
    
    Kept private to avoid cluttering global namespace.
    
    \param file_name_with_path_wrap Wrap for file_name_with_path
    \param log_file_wrap Wrap for log_file
    \param exe_name_wrap Wrap for exe_name
    \param user_name_wrap Wrap for user_name
    */
		inline void initialize_class( std::string &file_name_with_path_wrap, 
			std::ofstream &log_file_wrap, std::string &exe_name_wrap, std::string &user_name_wrap );
    
    /**
	  \brief The function used to do the actual writing onto the file
	
	  Kept private to avoid cluttering global namespace.

	  \param FreeText_wrap Wrap for FreeText
	  \param log_file_wrap Wrap for log_file
	  \param timer_wrap Wrap for timer
	  \param exe_name_wrap Wrap for exe_name
	  \param user_name_wrap Wrap for user_name
	  */
	  inline void writing_function( const std::string &FreeText_wrap, std::ofstream &log_file_wrap, 
			const std::string &exe_name_wrap, const std::string &user_name_wrap, bool isError = false );

	private:
    //! The file handler class
		std::ofstream log_file;
    //! The name of the executable calling the class
		std::string exe_name; 
    //! The current active user name
		std::string user_name;
    //! The free text to be written the log file. It is taken as input from user_name
		//std::string FreeText; 
    //! File path
    std::string file_name_with_path;
    //! Flag to initialize local logging
    bool consoleLogging;
    //! Flag to denote logging in GMT
    bool GMTLogging;
    //! This writes out user and exe name(s) to the log
    bool m_multiUserLog = false;
	};
}