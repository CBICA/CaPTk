#ifndef _UTILS_H
#define _UTILS_H

#include <fstream>
#include <iostream>
#include <iterator>
#include <ctime>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <itkNiftiImageIO.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>
#include <sstream>
#include <sys/stat.h>
using namespace std;
typedef int(*TXT_CALLBACK_TYPE)(const char*);

typedef float PixelType;
typedef itk::Image<PixelType, 3> ImageType;

namespace Utils
{

  inline string getBasePath(const std::string& fname)
  {
    size_t pos = fname.find_last_of("\\/");
    if (std::string::npos == pos)
    {
      return string();
    }
    else
    {
      return fname.substr(0, pos);
    }
  }
  inline bool dirExists(const std::string& path)
  {
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
    {
      return false;
    }
    else if (info.st_mode & S_IFDIR)
    {
      return true;
    }
    else
    {
      return false;
    }

  }
	inline bool fileExist(const std::string& name) {
		ifstream f(name.c_str());
		return f.good();
	}
	inline bool checkExtn(const std::string& name, const std::string& extn) {
		return  (name.substr(name.find_last_of(".") + 1) == extn);
	}
	static TXT_CALLBACK_TYPE g_listner = NULL;
	static string toString(const int& value)
	{
		std::ostringstream ostr;
		ostr << value;
		return  ostr.str(); 
	}

	template <typename T>
	static void writeCsv(const vector<T>& data, const string& fileName)
	{
		ofstream fs(fileName.c_str());
		if (fs.is_open())
		{
			for (size_t r = 0; r < data.size(); r++)
			{
				fs << data[r] << endl;
			}
			fs.close();
		}
		return;
	}
	template <typename T>
	static void writeCsv(const vector< vector<T> >& data, const string& fileName)
	{
		ofstream fs(fileName.c_str());
		if (fs.is_open())
		{
			for (size_t r = 0; r < data.size(); r++)
			{
				for (size_t c = 0; c < data[r].size(); c++)
				{
          fs << data[r][c];
          if (c != data[r].size() - 1)
          {
            fs << ", ";
          }
				}
				fs << endl;
			}
			fs.close();
		}
		return;
	}
	template <typename T>
	static  vector<vector< T > > readCsv(const string& filename)
	{
		vector< vector<T> > valuesOut;
		ifstream fin(filename.c_str());
		if (fin.is_open())
		{
			for (string line; getline(fin, line);)
			{
				replace(line.begin(), line.end(), ',', ' ');
				istringstream in(line);
				valuesOut.push_back(vector<T>(istream_iterator<T>(in), istream_iterator<T>()));
			}
		}
		return  valuesOut;
	}
	template <typename T>
	static std::vector<T> flatten(const std::vector< std::vector<T> >& v) {
		std::size_t total_size = 0;
		for (size_t index = 0; index != v.size(); ++index) {
			total_size += v[index].size();
		}
		std::vector<T> result;
		result.reserve(total_size);
		for (size_t index = 0; index != v.size(); ++index) {
			result.insert(result.end(), v[index].begin(), v[index].end());
		}
		return result;
	}
	template <typename T>
	static cv::Mat toMat(const vector<vector<T> > vecIn) {
		cv::Mat_<T> matOut(vecIn.size(), vecIn.at(0).size());
		for (int i = 0; i < matOut.rows; ++i) {
			for (int j = 0; j < matOut.cols; ++j) {
				matOut(i, j) = vecIn.at(i).at(j);
			}
		}
		return matOut;
	}
  template <typename T>
  static vector<vector<T> > toVector(cv::Mat matIn) 
  {
      vector<vector<T> > vecOut;
      for (int i = 0; i < matIn.rows; ++i) 
      {
        vector<T> row(matIn.ptr<T>(i), matIn.ptr<T>(i)+matIn.cols);
        vecOut.push_back(row);
      }
      return vecOut;
  }

	static std::clock_t    g_lastClock;
	static void resetTimer()
	{
		g_lastClock = std::clock();
	}
	static void setListner(const TXT_CALLBACK_TYPE listner)
	{
		g_listner = listner;
	}
	static void logText(const string msg = "")
	{
		if (g_listner != NULL)
		{
			g_listner(msg.c_str());
		}
		else
		{
			double sec = (std::clock() - g_lastClock) / (double)(CLOCKS_PER_SEC);
			stringstream ss;
			ss << msg << " Time: " << sec << " Sec (" << sec << "ms)" << std::endl;
			std::cout << ss.str();
			g_lastClock = std::clock();
		}
		return;
	}
	static void progressUpdate(const size_t iter, const size_t iterMax, const string& name)
	{
		if (iterMax == 0)
		{
			stringstream ss;
			ss << name << ":Progress= " << -1.0 << "%";
			logText(ss.str());
		}
		else if (iterMax < 100 || iter % (iterMax / 100) == 0)
		{
			stringstream ss;
			ss << name << ":Progress= " << (100.0*iter) / iterMax << "%";
			logText(ss.str());
		}
		return;
	}
	static ImageType::Pointer readNifti(const string& fileName)
	{
		typedef itk::ImageFileReader<ImageType> ReaderType;
		ImageType::Pointer image = ImageType::New();
		try
		{
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(fileName);
			reader->Update();
			image->Graft(reader->GetOutput());
		}
		catch (...)
		{
			Utils::logText("Exception in reading Nifti:" + fileName);
		}
		return image;
	}
	static void writeNifti(const string& inputFile, const string& outFile, ImageType::Pointer image)
	{
		itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFile.c_str(), itk::ImageIOFactory::ReadMode);
		imageIO->SetFileName(inputFile);
		imageIO->ReadImageInformation();
		itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
		nifti_io->SetPixelType(imageIO->GetPixelType());
		//TBD
	}
	template <class T>
	static T unpack(uint8_t* packed, bool littleEndian = true) {
		if (littleEndian) {
			return *reinterpret_cast<T*>(packed);
		}
		else {
			uint8_t backwards[sizeof(T)];
			for (uint8_t i = 0; i < sizeof(T); i++) {
				backwards[i] = packed[sizeof(T) - i - 1];
			}
			return *reinterpret_cast<T*>(backwards);
		}
	}
	template <class T>
	static void pack(T* unpacked, uint8_t* packed, bool littleEndian = true) {
		if (littleEndian) {
			for (uint8_t i = 0; i < sizeof(T); i++) {
				packed[i] = reinterpret_cast<uint8_t*>(unpacked)[i];
			}
		}
		else {
			for (uint8_t i = 0; i < sizeof(T); i++) {
				packed[i] = reinterpret_cast<uint8_t*>(unpacked)[sizeof(T) - i - 1];
			}
		}
	}
	static float readVal(ifstream& f)
	{
		uint8_t byte[4];
		f.read((char *)byte, 4);
		return unpack<float>(byte, false);
	}
	static void writeVal(ofstream& f, float val)
	{
		uint8_t byte[4];

		pack(&val, byte, false);
		f.write((char *)byte, 4);
		return;
	}

	static vector<cv::Point3f> readFiberBfloat(const string& bfloatFile, vector<long>& idxs)
	{

		idxs.clear();
		vector<cv::Point3f> points;
		ifstream fs(bfloatFile.c_str(), ios::binary);
		if (fs.is_open())
		{
			fs.seekg(0, std::ios_base::end);
			size_t sizeMax = fs.tellg();
			fs.seekg(0, std::ios_base::beg);
			size_t size = 0;
			while (size<sizeMax)
			{

				int numPoints = Utils::readVal(fs);
				int seedPoint = Utils::readVal(fs);
        seedPoint = seedPoint;
				assert(numPoints >= 0);
				size = size + (3 * numPoints + 2)* sizeof(float);
				assert(size <= sizeMax);
				for (size_t i = 0; i < (uint)numPoints; i++)
				{
					cv::Point3f pt;
					pt.x = Utils::readVal(fs);
					pt.y = Utils::readVal(fs);
					pt.z = Utils::readVal(fs);
					points.push_back(pt);
				}
				idxs.push_back(numPoints);
			}
			fs.close();
		}
		return points;
	}

	static void threshTruncVal(cv::Mat& data, const double thresh, const double val)
	{
		cv::MatIterator_<double> it, end;
		for (it = data.begin<double>(), end = data.end<double>(); it != end; ++it)
		{
			if ((*it) <= thresh)
			{
				(*it) = val;
			}
		}
		return;
	}
	static void dfIdfTransform(cv::Mat& X)
	{
		int n_samples = X.rows;
		vector<double> df(X.cols);
		for (int c = 0; c < X.cols; c++)
		{
			df[c] = countNonZero(X.col(c)) + 1;//+1 to avoid div by zero
		}
		n_samples += 1;
		for (int i = 0; i < X.cols; i++)
		{
			df[i] = log(double(n_samples) / df[i]) + 1.0; //Inverse of idf
		}
		cv::Mat idf_diag = cv::Mat::diag(cv::Mat(df));
		X = X* idf_diag;
		return;
	}
	static double gammaLn(const double xx)
	{
		double x, y, tmp, ser;
		static double cof[6] = { 76.18009172947146, -86.50532032941677,
			24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
			-0.5395239384953e-5 };
		y = xx;
		x = xx;
		tmp = x + 5.5;
		tmp -= (x + 0.5)*log(tmp);
		ser = 1.000000000190015;
		for (int j = 0; j <= 5; j++)
		{
			ser += cof[j] / ++y;
		}
		return(-tmp + log(2.5066282746310005 * ser / x));

	}
	static cv::Mat gammaLn(const cv::Mat _data)
	{
		cv::Mat data = _data.clone();
		for (int i = 0; i < data.rows; i++)
		{
			for (int j = 0; j < data.cols; j++)
			{
				data.at<double>(i, j) = gammaLn(data.at<double>(i, j));
			}
		}
		return data;
	}
}
#endif
