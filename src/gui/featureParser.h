#ifndef _featureParser_h_
#define _featureParser_h_

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#if defined(_WIN32)&& !defined(__GNUC__)
  #pragma warning( disable : 4503 )// Warnig related to decorated name length size
#endif
using namespace std;
//TBD add documentation once class is frozen
inline vector<string> splitStr(const string& str, const string& delim)
{
  vector<string> tokens;
  size_t prev = 0, pos = 0;
  do
  {
    pos = str.find(delim, prev);
    if (pos == string::npos) pos = str.length();
    string token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}
inline void removeChar(string& str, char ch)
{
  str.erase(remove(str.begin(), str.end(), ch), str.end());
}
class FeatureParser
{
public:
  FeatureParser(const string& fileName)
  {
    m_fileName = fileName;
    m_featureMap = prepareFeatureMap();
  }
  ~FeatureParser()
  {
  }
  string getLastError()
  {
    return m_lastError;
  }
  map<string, vector<map<string, string> > > getFeatureMap()
  {
    return m_featureMap;
  }
  void setFeatureMap(map<string, vector<map<string, string> > > featureMap)
  {
    m_featureMap = featureMap;
  }
private:
  map<string, vector<map<string, string> > > prepareFeatureMap()
  {
    vector<string> header;
    vector<vector<string> > matrix;
    readCsv(m_fileName, header, matrix);
    map < string, vector<map<string, string> > > features;
    for (auto row : matrix)
    {
      map<string, string> props;
      for (size_t i = 1; i < row.size(); i++)
      {
        props[header[i]] = row[i];
        if (header[i] == "Default")
        {
          props["Value"] = row[i];
        }
      }
      features[row[0]].push_back(props);
    }
    return features;
  }
  string trim(const string &s)
  {
    auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c){return std::isspace(c); });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c){return std::isspace(c); }).base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
  }
  bool readCsv(const string& fileName, vector<string>& header, vector<vector<string> >&   matrix)
  {
    std::ifstream       file(fileName.c_str());
    if (!file)
    {
      m_lastError = "Could not read: " + fileName;
      return false;
    }
    while (file)
    {
      string line;
      getline(file, line);
      stringstream lineStream(line);
      vector<string> row;
      string cell;
      while (getline(lineStream, cell, ','))
      {
        row.push_back(trim(cell));
      }
      if (!row.empty() && header.empty())
      {
        header = row;
      }
      else if(!row.empty())
      {
        row.resize(header.size());
        if (row[0].empty() && !matrix.empty())
        {
          row[0] = matrix.back()[0]; //Copy previous feature family name if empty
        }
        matrix.push_back(row);
      }
    }
    if (header.empty() || matrix.empty())
    {
      m_lastError = "Empty file: " + fileName;
      return false;
    }
    return true;
  }

  string m_lastError;
  string m_fileName;
  map<string, vector<map<string, string> > > m_featureMap;
};


#endif