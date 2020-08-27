#include "testWhiteStripe.h"
#include <ctime>


//------------------Utility functions----------------------------------
inline bool icasecmp(const string& l, const string& r)
{
  return l.size() == r.size()
    && equal(l.cbegin(), l.cend(), r.cbegin(),
    [](string::value_type l1, string::value_type r1)
  { return toupper(l1) == toupper(r1); });
}
inline bool cmdOptionExists(const vector<string>& allArgs, const string& option)
{
  for (auto &arg : allArgs)
  {
    if (icasecmp(arg, option))
    {
      return true;
    }
  }
  return false;
}
inline string getCmdOption(const vector<string>& allArgs, const std::string & option)
{
  for (size_t i = 0; i < allArgs.size() - 1; i++)
  {
    if (icasecmp(allArgs[i], option))
    {
      return allArgs[i + 1];
    }
  }
  return string();
}
inline bool fileExist(const std::string& name) {
  ifstream f(name.c_str());
  return f.good();
}
void showVersionInfo()
{
  cout << "WhiteStripe version: " << WHITESTRIPE_VERSION << endl;
}
inline void copyrightNotice()
{
  cout <<
    "\n==========================================================================\n" <<
    "Contact: software@cbica.upenn.edu\n\n" <<
    "Copyright (c) 2017 University of Pennsylvania. All rights reserved.\n" <<
    "See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html" <<
    "\n==========================================================================\n";
}
void showhelpinfo(const string exeName)
{
  cout << "Usage:   " << exeName << " [-option] [argument]" << endl;
	cout << "Option:  " << "-h show help information" << endl;
  cout << "         " << "-v show version information" << endl;
	cout << "         " << "-i inputFile" << endl;
	cout << "         " << "-o outputFile[optional]" << endl;
	cout << "         " << "-t fileType: 0(T1) 1(T2)[optional, default 0]" << endl;
	cout << "         " << "-m maxTissues: 5 or 10 or 20[optional, default 5]" << endl;
  cout << "         " << "-zStart start z slice number for cropping[optional, default 80]" << endl;
  cout << "         " << "-zStop stop z slice number for cropping[optional, default 120]" << endl;
  cout << "Example: " << exeName << " -i inFile.nii.gz -t 0" << endl;
}
int runWhiteStripeCmd(const vector<string>& allArgs)
{
	std::clock_t    start;
	start = std::clock();
	try
	{
		string ansStr, inFile, outFile= "out.nii.gz";
		float twsWidth = 0.05;
		int sliceStartZ = 80;
		int sliceStopZ = 120;
		int tissuesMax = 5;
		float smoothMax = 10.0;
		float smoothDelta = 0.5;
		int histSize = 2000;
		bool bT1 = true;

    
    ansStr = getCmdOption(allArgs, "-i");
		if (!ansStr.empty())
		{
			inFile = ansStr;
		}
		else
		{
			cout << "Error.No input file specified:- "<< endl;
      showhelpinfo(allArgs[0]);
      return EXIT_FAILURE;
		}
    ansStr = getCmdOption(allArgs, "-o");
		if (!ansStr.empty())
		{
			outFile = ansStr;
		}

    ansStr = getCmdOption(allArgs, "-t");
		if (!ansStr.empty())
		{
			bT1 = (0 == atoi(ansStr.c_str()));
		}

    ansStr = getCmdOption(allArgs, "-m");
		if (!ansStr.empty())
		{
			tissuesMax =atoi(ansStr.c_str());
		}

    ansStr = getCmdOption(allArgs, "-zStart");
    if (!ansStr.empty())
    {
      sliceStartZ = atoi(ansStr.c_str());
    }
    ansStr = getCmdOption(allArgs, "-zStop");
    if (!ansStr.empty())
    {
      sliceStopZ= atoi(ansStr.c_str());
    }

		ImageType::Pointer img, mask;
		img = Utils::readNifti(inFile);
		if (img.IsNull())
		{
			cout << "Error reading:- "<< inFile << endl;
      return EXIT_FAILURE;
		}


		WhiteStripe obj;
		obj.setParams(twsWidth, sliceStartZ, sliceStopZ, tissuesMax, smoothMax, smoothDelta, histSize, bT1);

		ImageType::Pointer normImg = obj.process(img, mask);
		Utils::writeNifti(inFile, outFile, normImg);
		cout << "Processed-" << inFile << " in: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
	}
	catch (...)
	{
		cout << "Exception in processing  " << endl;
    return EXIT_FAILURE;
	}

  return EXIT_SUCCESS;
}
int runWhiteStripeTest(const vector<string>& allArgs)
{

  try
  {
    if (cmdOptionExists(allArgs, "Quantile"))
    {
      return TestWhiteStripe().testQuantile();

    }
    if (cmdOptionExists(allArgs, "ImgMean"))
    {
      return TestWhiteStripe().testImgMean();
    }
    if (cmdOptionExists(allArgs, "ImgMean"))
    {
      return TestWhiteStripe().testWhiteStripe();
    }

  }
  catch (...)
  {
    cout << "Exception in processing  " << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
int main(int argc, char *argv[]) {
  copyrightNotice();//TBD remove this if not needed 
  std::vector<string> allArgs(argv, argv + argc);
  if (cmdOptionExists(allArgs, "-v") || cmdOptionExists(allArgs, "-version"))
  {
    showVersionInfo();
    return EXIT_SUCCESS;
  }
  if (cmdOptionExists(allArgs, "-h") || cmdOptionExists(allArgs, "-help"))
  {
    showhelpinfo(allArgs[0]);
    return EXIT_SUCCESS;
  }
  if (cmdOptionExists(allArgs, "-test"))
  {
    return runWhiteStripeTest(allArgs);
  }
  else
  {
    return runWhiteStripeCmd(allArgs);
  }
  return EXIT_SUCCESS;
}
