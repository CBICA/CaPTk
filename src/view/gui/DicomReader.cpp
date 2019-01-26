#include "DicomReader.h"

DicomReader::DicomReader()
{
  reader = ReaderType::New();
  dicomIO = ImageIOType::New();
  nameGenerator = NamesGeneratorType::New();
  reader->SetImageIO(dicomIO);
  tagvalueMap = MapContainerType::New();
}

DicomReader::~DicomReader()
{
}

void DicomReader::SetDirectoryPath(std::string path)
{
  nameGenerator->SetInputDirectory(path);
}

void DicomReader::ReadMetaData()
{
  tagvalueMap->Initialize(); //clear the map
  fileNames = nameGenerator->GetInputFileNames();
  reader->SetFileNames(fileNames);

  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cout << ex << std::endl;
  }

  dictionary = dicomIO->GetMetaDataDictionary();

  DictionaryType::ConstIterator itr = dictionary.Begin();
  DictionaryType::ConstIterator end = dictionary.End();

  while (itr != end)
  {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;

    MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>(entry.GetPointer());

    std::string label;
    std::pair<std::string, std::string> labelValuePair;
    if (entryvalue)
    {
      std::string tagkey = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      bool found = itk::GDCMImageIO::GetLabelFromTag(tagkey, label);
      if (found)
      {
        labelValuePair.first = label;
        labelValuePair.second = tagvalue;
        tagvalueMap->InsertElement(tagkey, labelValuePair);
      }

      //std::cout << tagkey <<  " = " << label << " " << tagvalue << std::endl;
    }

    ++itr;
  }

}

void DicomReader::PrintMetaData()
{
  MapContainerType::Iterator itr;
  std::pair<std::string, std::string> labelValuePair;
  for (itr = tagvalueMap->Begin(); itr != tagvalueMap->End(); ++itr)
  {
    labelValuePair = itr->Value();
    std::cout << itr->Index() << " " << labelValuePair.first << " " << labelValuePair.second << std::endl;
  }

}

std::map<std::string, std::pair<std::string, std::string>> DicomReader::GetMetaDataMap()
{
  return tagvalueMap->CastToSTLContainer();

}

bool DicomReader::GetTagValue(std::string tag, std::string &label, std::string &value)
{
  bool status = false;
  DictionaryType::ConstIterator tagItr = dictionary.Find(tag);
  if (tagItr == dictionary.End())
    status = false;
  else
  {
    MetaDataStringType::ConstPointer entryvalue = dynamic_cast<const MetaDataStringType *>(tagItr->second.GetPointer());
    if (entryvalue)
    {
      //status = true;
      value = entryvalue->GetMetaDataObjectValue();
      status = itk::GDCMImageIO::GetLabelFromTag(tag, label);
    }
    else
    {
      status = false; //entry was not string
    }
  }

  return status;
}
