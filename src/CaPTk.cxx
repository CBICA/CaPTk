/**
\file CAPTk.cpp

\brief Main entry point for CAPTk
*/

#include <QApplication>
#include "fMainWindow.h"
//#include "CAPTk.h"

#include "vtkFileOutputWindow.h"
#include "itkFileOutputWindow.h"
//#include "vtkOutputWindow"

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "yaml-cpp/yaml.h"

#include "CheckOpenGLVersion.h"

///// debug
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
//#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )      
//#define new DBG_NEW   
///// debug

void echoCWLFiles(std::vector< std::string > inputCWLFiles)
{
  std::cout << "Availble CWL applications (refer to individual CLI usage): \n\n";
  for (size_t i = 0; i < inputCWLFiles.size(); i++)
  {
    auto cwlFileBase = cbica::getFilenameBase(inputCWLFiles[i]);
    std::cout << "\t" << cwlFileBase << "\n";
  }
}

#ifdef _WIN32
// ensures no console pops up when launching the program
int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE prevInstance, LPSTR lpCmdLine, int nShowCmd)
{
  int argc = __argc;
  char **argv = __argv;
#else
int main(int argc, char** argv)
{
#endif

  std::string cmd_inputs, cmd_mask, cmd_tumor, cmd_tissue;
  float cmd_maskOpacity = 1;
  bool comparisonMode = false;

  // this is used to populate the available CWL files for the cli  
  auto cwlFolderPath = 
#ifdef CAPTK_PACKAGE_PROJECT
  cbica::normPath(cbica::getExecutablePath() + 
#ifdef __APPLE__
    "../Resources/etc/cwlDefinitions/"
#else
    "../etc/cwlDefinitions/"
#endif
    )
#else
  std::string(PROJECT_SOURCE_DIR) + "/data/cwlFiles"
#endif
  ;
  // std::cout << cwlFolderPath + "\n"; 

  auto cwlFiles = cbica::filesInDirectory(cwlFolderPath);

  // parse the command line
  auto parser = cbica::CmdParser(argc, argv, "CaPTk");
  parser.ignoreArgc1();

  parser.addOptionalParameter("i", "images", cbica::Parameter::FILE, "NIfTI or DICOM", "Input coregistered image(s) to load into CaPTk", "Multiple images are delineated using ','");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "NIfTI or DICOM", "Input mask [coregistered with image(s)] to load into CaPTk", "Accepts only one file");
  parser.addOptionalParameter("mo", "maskOpacity", cbica::Parameter::FLOAT, "0-1", "Opacity of the Input mask", "Needs 'm' to be passed");
  parser.addOptionalParameter("tu", "tumorPt", cbica::Parameter::FILE, ".txt", "Tumor Point file for the image(s) being loaded");
  parser.addOptionalParameter("ts", "tissuePt", cbica::Parameter::FILE, ".txt", "Tissue Point file for the image(s) being loaded");
  parser.addOptionalParameter("a", "advanced", cbica::Parameter::BOOLEAN, "none", "Advanced visualizer which does *not* consider", "origin information during loading");
  parser.addOptionalParameter("c", "comparisonMode", cbica::Parameter::BOOLEAN, "true or false", "Enable/Disable comparison mode", "comparison mode during loading");

  //parser.exampleUsage("-i C:/data/input1.nii.gz,C:/data/input2.nii.gz -m C:/data/inputMask.nii.gz -tu C:/data/init_seed.txt -ts C:/data/init_GLISTR.txt");
  parser.addExampleUsage("-i C:/data/input1.nii.gz,C:/data/input2.nii.gz -m C:/data/inputMask.nii.gz -tu C:/data/init_seed.txt -ts C:/data/init_GLISTR.txt",
    "Load the input images and ROI with seed points");
  parser.addApplicationDescription("Entry point for all CaPTk applications");

  // check for CWL command coming in through the command line after "CaPTk"
  if (argc > 1)
  {
    for (auto & file : cwlFiles)
    {
      auto cwlFileBase = cbica::getFilenameBase(file);
      auto cwlFileBase_actual = cwlFileBase;
      std::transform(cwlFileBase.begin(), cwlFileBase.end(), cwlFileBase.begin(), ::tolower);
      auto argv_1 = std::string(argv[1]);
      argv_1 = cbica::getFilenameBase(argv_1, false);
      std::transform(argv_1.begin(), argv_1.end(), argv_1.begin(), ::tolower);

      // Check for filename without cwl extension
      if (cwlFileBase.find(argv_1) != std::string::npos)
      {
        // Get base command
        //std::ofstream selected_file;
        //selected_file.open(file.c_str());
        auto config = YAML::LoadFile(file);
        // Get all args passed to application
        std::string argv_complete;
        for (size_t i = 2; i < argc; i++) // 2 because the argv[1] is always the "application"
        {
          argv_complete += " \"" + std::string(argv[i]) + "\""; // add double quote for command tokenization
        }
        // Pass them in
        auto commandToRun = getApplicationPath(cwlFileBase_actual) + argv_complete;
        //std::cout << "[DEBUG] commandToRun: " << commandToRun << "\n";
// #ifndef WIN32 
        return std::system(commandToRun.c_str());
// #else
//         auto returnCode = std::system(commandToRun.c_str());
//         std::system("pause");
//         return returnCode;
// #endif
      }
    }
  }

  if (parser.isPresent("u", false))
  {
    parser.echoUsage();
    echoCWLFiles(cwlFiles);
    return EXIT_SUCCESS;
  }

  if (parser.isPresent("h", false))
  {
    parser.echoHelp();
    echoCWLFiles(cwlFiles);
    return EXIT_SUCCESS;
  }

  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", cmd_inputs);
  }
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", cmd_mask);
    if (parser.isPresent("mo"))
    {
      parser.getParameterValue("m0", cmd_maskOpacity);
    }
  }
  if (parser.isPresent("tu"))
  {
    parser.getParameterValue("tu", cmd_tumor);
  }
  if (parser.isPresent("ts"))
  {
    parser.getParameterValue("ts", cmd_tissue);
  }
  if (parser.isPresent("c"))
  {
    parser.getParameterValue("c", comparisonMode);
  }

#if defined(__linux__)
  //auto defaultFormat = QVTKOpenGLWidget::defaultFormat();
  //// defaultFormat.setSamples(0);
  //defaultFormat.setVersion(3, 0);
  //QSurfaceFormat::setDefaultFormat(defaultFormat);
#else
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
#endif

  // high DPI fixes
  QApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
  QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
  //QCoreApplication::setAttribute(Qt::AA_UseOpenGLES);

  QApplication app(argc, argv);

  //cbica::setEnvironmentVariable("QT_QPA_PLATFORM_PLUGIN_PATH", captk_currentApplicationPath + "/platforms");
  //cbica::setEnvironmentVariable("QT_OPENGL", "software");


  // starting the OpenGL version checking 
  const std::string openGLVersionCheckFile = loggerFolderBase + "openglVersionCheck.txt";
  if (!cbica::isFile(openGLVersionCheckFile))
  {
    std::cout << "Checking for compatible OpenGL - this will happen only once.\n";
    std::string msg;
    bool minimumVersionNotFound = false;
#if WIN32
    CheckOpenGLVersion checker(hInstance);
#else
    CheckOpenGLVersion checker;
#endif
    if (!checker.hasVersion_3_2())
    {
      std::string msg = "A working 3.2 version of OpenGL was not found in your hardware/software combination; consequently, CaPTk's GUI will not work; all CLIs will work as expected.\n\n";
      msg += "\tOpenGL Version : " + checker.version + "\n";
      msg += "\tOpenGL Renderer: " + checker.renderer + "\n";
      msg += "\tOpenGL Vendor  : " + checker.vendor + "\n\n";
	  msg += "Please install OpenGL version 3.2 or greater, or use the command line interface to run CaPTk.";
	  msg += "\n\nCheck the documentation for details.";
	  ShowErrorMessage(msg);
#if WIN32
      cbica::sleep(1000);
      return EXIT_FAILURE;
#else // Attempt software rendering on non-Windows platforms, warn user 
      cbica::setEnvironmentVariable("QT_OPENGL", "software");
	  std::string softwareRenderingMsg = "WARNING: Trying to run CaPTk GUI using software rendering - this might not work on all systems and in those cases, only the CLI will be available.\n";
	  std::cerr << softwareRenderingMsg;
	  ShowErrorMessage(softwareRenderingMsg);
#endif
    }
    else
    {
      std::cout << "Compatible OpenGL was found. This check will not happen again for this machine.\n";
      std::ofstream myFile;
      myFile.open(openGLVersionCheckFile.c_str());
      myFile << "Compatible OpenGL version present.\n";
      myFile.close();
    }
  }

  //! redirect the vtkoutputwindow contents to file
  auto fileOutputWindow = vtkSmartPointer< vtkFileOutputWindow >::New();
  fileOutputWindow->SetFileName((loggerFolderBase + "vtk_errors.txt").c_str());
  vtkOutputWindow::SetInstance(fileOutputWindow);

  //! redirect the itk output window contents to file
  auto itkOutputWindow = itk::FileOutputWindow::New();
  itkOutputWindow->SetFileName((loggerFolderBase + "itk_errors.txt").c_str());
  itkOutputWindow->FlushOn();
  itk::OutputWindow::SetInstance(itkOutputWindow);

  fMainWindow window; // initialize main app
  if (parser.isPresent("a"))
  {
    window.EnableAdvancedVisualizer();
  }

  // get everything to QString
  std::vector< QString > inputFiles_QString;
  QString inputMask;

  if (!cmd_inputs.empty()) // the check is only for input images since images can be rendered with ROI but not the other way round
  {
    auto inputFiles = cbica::stringSplit(cmd_inputs, ",");
    auto inputMasks = cbica::stringSplit(cmd_mask, ",");

    for (size_t i = 0; i < inputFiles.size(); i++)
    {
      if (cbica::fileExists(inputFiles[i])) // check for file validity
      {
        inputFiles_QString.push_back(inputFiles[i].c_str());
      }
    }

    if (inputFiles_QString.size() != 0)
    {
      for (size_t i = 0; i < inputMasks.size(); i++)
      {
        if (cbica::fileExists(inputMasks[i])) // check for file validity
        {
          inputMask = inputMasks[i].c_str();
          break; // because only 1 annotated ROI can be drawn up at one time
        }
      }
    }
    else // at this point, give an error if there are no valid image files passed
    {
      std::cerr << "No valid images were provided.\n";
      return EXIT_FAILURE;
    }
  }

  cbica::createDir(loggerFolderBase);
  cbica::createDir(loggerFolder);
  //cbica::createDir(captk_StuffFolderBase);
  //cbica::createDir(captk_SampleDataFolder);
  //cbica::createDir(captk_PretrainedFolder);

#ifndef _WIN32
  std::string old_locale = setlocale(LC_NUMERIC, NULL);
  setlocale(LC_NUMERIC, "POSIX");
#endif

  app.processEvents();

  if (inputFiles_QString.size() == 0)
  {
#ifndef CAPTK_PACKAGE_PROJECT
    std::string testData = "/data/AAAC0_flair_pp_shrunk.nii.gz";
    std::string testData_mask = "/data/AAAC0_flair_pp_shrunk_testTumor.nii.gz";

    // if test data is present, load it up otherwise forget about it
    std::string appDir = captk_currentApplicationPath;
    if (cbica::fileExists(appDir + testData))
    {
      inputFiles_QString.push_back((appDir + testData).c_str());
      inputMask.push_back((appDir + testData_mask).c_str());
    }
    else if (cbica::fileExists(appDir + "/.." + testData))
    {
      inputFiles_QString.push_back((appDir + "/.." + testData).c_str());
      inputMask.push_back((appDir + "/.." + testData_mask).c_str());
    }
    else if (cbica::fileExists(appDir + "/../.." + testData))
    {
      inputFiles_QString.push_back((appDir + "/../.." + testData).c_str());
      inputMask.push_back((appDir + "/../.." + testData_mask).c_str());
    }
#endif
  }
  else
    window.loadFromCommandLine(inputFiles_QString, comparisonMode, inputMask.toStdString(),
      cmd_maskOpacity, cmd_tumor, cmd_tissue/*, true*/); // at this point, inputFiles_QString will have at least 1 value



  //window.setFixedSize(QSize(2735, 1538)); // useful when doing video recording from SP4 and maintain 16x9 ratio
  window.showMaximized();
  window.show();

  // show the "about" screen in the first run
  if (!cbica::fileExists(tutorialScreen))
  {
    // auto rec = QApplication::desktop()->screenGeometry();
    // std::cout << "Detected Size: " << rec.width() << "x" << rec.height() << "\n";
    window.about();
  }



#ifndef _WIN32
  setlocale(LC_NUMERIC, old_locale.c_str());
#endif

  return app.exec();
}

