/**
\file  cbicaLogging.cpp

\brief Implementation of the Logging class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2016 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#include "cbicaLogging.h"

namespace cbica
{ 
  Logging::Logging(const std::string file_name, const std::string FreeText_input = "") 
  {
    consoleLogging = false; // console output is disabled if file is provided
    GMTLogging = false;
    file_name_with_path = file_name;
    initialize_class(file_name_with_path, log_file, exe_name, user_name);
    if (FreeText_input.empty())
    {
      writing_function("", log_file, exe_name, user_name);
    }
    else
    {
      writing_function("" + FreeText_input, log_file, exe_name, user_name);
    }
    log_file.close();
  }
  
  Logging::Logging() 
  {	
    file_name_with_path = "";
    consoleLogging = true; // no file path translates to console output
    GMTLogging = false;
    initialize_class(file_name_with_path, log_file, exe_name, user_name); 
    //writing_function("", log_file, exe_name, user_name); 
  }

  void Logging::UseNewFile(const std::string &newLogFile = "")
  {
    consoleLogging = false;
    GMTLogging = false;
    if (newLogFile.empty())
    {
      file_name_with_path = cbica::createTmpDir() + cbica::getExecutableName() + "-log.txt";
    }
    else
    {
      file_name_with_path = newLogFile;
    }
    initialize_class(file_name_with_path, log_file, exe_name, user_name);
    //writing_function("", log_file, timer, exe_name, user_name);
    log_file.close();
  }

  Logging::Logging(const Logging &origin)
  {
    //file_name_with_path = origin.file_name_with_path;
  }

  Logging::~Logging() 
  { 
    //log_file.close(); // this is no longer needed as the Write() command opens and closes the log_file every single time
    // this is done so that the OS don't keep a file handle lock on the log file and the user can visualize the log as it happens
  }

  void Logging::EnableTextLogging(const std::string &newLogFile = "")
  {
    UseNewFile(newLogFile);
  }

  void Logging::EnableGMTLogging()
  {
    GMTLogging = true;
  }

  void Logging::EnableConsoleLogging()
  {
    consoleLogging = true;
    initialize_class(file_name_with_path, log_file, exe_name, user_name);
    //writing_function("", log_file, timer, exe_name, user_name);
    log_file.close();
  }

  void Logging::WriteError(const std::string FreeText_input)
  {
    // assumes file exists because constructor writes the file once
    log_file.open(file_name_with_path.c_str(), std::ios_base::app);
    writing_function("" + FreeText_input, log_file, exe_name, user_name, true);
    log_file.close();
  }

  void Logging::Write(const std::string FreeText_input = "")
  { 
    // assumes file exists because constructor writes the file once
    if (!file_name_with_path.empty())
    {
      log_file.open(file_name_with_path.c_str(), std::ios_base::app);
      while (!log_file.is_open())
      {
        cbica::sleep(100);
        log_file.open(file_name_with_path.c_str(), std::ios_base::app);
      }
    }

    if (FreeText_input.empty())
    {
      writing_function("", log_file, exe_name, user_name);
    }
    else
    {
      writing_function("" + FreeText_input, log_file, exe_name, user_name);
    }
    log_file.close();
  }
  
	inline void Logging::initialize_class(std::string &file_name_with_path_wrap, std::ofstream &log_file_wrap, 
		std::string &exe_name_wrap, std::string &user_name_wrap )
	{
    if (!consoleLogging)
    {
      if (file_name_with_path_wrap == "")
      {
        file_name_with_path_wrap = cbica::createTmpDir() + "temp.log";
        std::cout << "Blank file name provided. A new file has been created with path: "
          << file_name_with_path_wrap << std::endl;
      }
      if (cbica::fileExists(file_name_with_path_wrap))
      {
        //std::cout << "File name '" << file_name_with_path_wrap << "' has been found. Appending.\n";
        log_file_wrap.open(file_name_with_path_wrap.c_str(), std::ios_base::app); // append to existing file
      }
      else
      {
        log_file_wrap.open(file_name_with_path_wrap.c_str());
      }
    }

    if (exe_name_wrap == "")
    {
      exe_name_wrap = cbica::getExecutableName();
    }
    if (user_name_wrap == "")
    {
      user_name_wrap = cbica::getUserName();
    }
	}

	inline void Logging::writing_function( const std::string &FreeText_wrap, std::ofstream &log_file_wrap, 
		const std::string &exe_name_wrap, const std::string &user_name_wrap, bool isError )
  {
    // obtain current time and add exeName, userName along with proper delineation
    std::string timeExeUser;
    if (GMTLogging)
    {
      timeExeUser = cbica::getCurrentGMTDateAndTime() + ";" + exe_name_wrap + ";" + user_name_wrap + ";";
    }
    else
    {
      timeExeUser = cbica::getCurrentLocalDateAndTime() + ";" + exe_name_wrap + ";" + user_name_wrap + ";";
    }

    if (consoleLogging)
    {
      if (isError)
      {
        std::cerr << timeExeUser + FreeText_wrap << "\n";
      }
      else
      {
        std::cout << timeExeUser + FreeText_wrap << "\n";
      }
    }
    else
    {
      log_file_wrap << timeExeUser.c_str() << FreeText_wrap.c_str() << "\n";
    }
	}

  std::string Logging::getLoggingFileName()
  {
    return file_name_with_path;
  }

} // end namespace cbica