#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaCmdParser.h"
#include "cbicaITKImageInfo.h"


int main(int argc, char **argv)
{
  std::cerr << "[DEPRECATED] This application can only be run via the GUI.\n";
  return EXIT_SUCCESS;
  //auto parser = cbica::CmdParser(argc, argv, "DirectionalityEstimate");
  //parser.addRequiredParameter("l", "labelMap", cbica::Parameter::FILE, "NIfTI file", "The label map on which calculations are to be done");
  //parser.addRequiredParameter("o", "outputFile", cbica::Parameter::FILE, "text file", "Where output gets stored");
  //parser.addOptionalParameter("r", "real", cbica::Parameter::FLOAT, "labelMap Dimensions", "The real world coordinates (in mm) of file",
  //  "Take example from tissue point file", "Needs to be in same dimensionality as labelMap", "Delineation is done using ','");
  //parser.addOptionalParameter("i", "index", cbica::Parameter::STRING, "labelMap Dimensions", "The image indeces (in voxels) of file",
  //  "Needs to be in same dimensionality as labelMap", "Delineation is done using ','");

  //std::string file_labelMap, file_output;

  //parser.getParameterValue("l", file_labelMap);
  //parser.getParameterValue("o", file_output);
  //parser.addApplicationDescription("This does a directionaliy estimation on the pre and post contrast injections for brains");
  //parser.addExampleUsage("-l c:/img_timepoint_1.nii.gz -o c:/output.txt", "This does directionality calculation based on the inputs");


  //std::string output_message = "Output Point:\nCoordinate = ";
  //std::string seedPoint;
  //bool realCoordinatesPassed = false;

  //if (parser.isPresent("r"))
  //{
  //  realCoordinatesPassed = true;
  //  parser.getParameterValue("r", seedPoint);
  //}

  //if (!realCoordinatesPassed)
  //{
  //  if (parser.isPresent("i"))
  //  {
  //    parser.getParameterValue("i", seedPoint);
  //  }
  //}

  //auto seedPoint_vec = cbica::stringSplit(seedPoint, ",");

  //auto imageInfo = cbica::ImageInfo(file_labelMap);

  //switch (imageInfo.GetImageDimensions())
  //{
  //case 2:
  //{
  //  using ImageType = itk::Image< float, 2 >;

  //  auto inputLabel = cbica::ReadImage< ImageType >(file_labelMap);

  //  ImageType::IndexType seedPoints_index, output_index;

  //  auto imageOrigin = inputLabel->GetOrigin();
  //  auto imageSpacing = inputLabel->GetSpacing();

  //  if (realCoordinatesPassed)
  //  {
  //    seedPoints_index[0] = std::abs((std::atof(seedPoint_vec[0].c_str()) - imageOrigin[0]) / imageSpacing[0]);
  //    seedPoints_index[1] = std::abs((std::atof(seedPoint_vec[1].c_str()) - imageOrigin[1]) / imageSpacing[1]);
  //  }
  //  else
  //  {
  //    seedPoints_index[0] = std::atof(seedPoint_vec[0].c_str());
  //    seedPoints_index[1] = std::atof(seedPoint_vec[1].c_str());
  //  }

  //  auto output = cbica::GetMaxDistanceInLabelMap< ImageType >(inputLabel, seedPoints_index);
  //  output_index = output.second;

  //  output_message += "[" + std::to_string(output_index[0]) + "," + std::to_string(output_index[1]) + "]"
  //    + "; Distance = " + std::to_string(output.first) + "\n";

  //  break;
  //}
  //case 3:
  //{
  //  using ImageType = itk::Image< float, 3 >;

  //  auto inputLabel = cbica::ReadImage< ImageType >(file_labelMap);

  //  ImageType::IndexType seedPoints_index, output_index;

  //  auto imageOrigin = inputLabel->GetOrigin();
  //  auto imageSpacing = inputLabel->GetSpacing();

  //  if (realCoordinatesPassed)
  //  {
  //    seedPoints_index[0] = std::abs((std::atof(seedPoint_vec[0].c_str()) - imageOrigin[0]) / imageSpacing[0]);
  //    seedPoints_index[1] = std::abs((std::atof(seedPoint_vec[1].c_str()) - imageOrigin[1]) / imageSpacing[1]);
  //    seedPoints_index[2] = std::abs((std::atof(seedPoint_vec[2].c_str()) - imageOrigin[2]) / imageSpacing[2]);
  //  }
  //  else
  //  {
  //    seedPoints_index[0] = std::atof(seedPoint_vec[0].c_str());
  //    seedPoints_index[1] = std::atof(seedPoint_vec[1].c_str());
  //    seedPoints_index[2] = std::atof(seedPoint_vec[2].c_str());
  //  }

  //  auto output = cbica::GetMaxDistanceInLabelMap< ImageType >(inputLabel, seedPoints_index);
  //  output_index = output.second;

  //  output_message += "[" + std::to_string(output_index[0]) + "," + std::to_string(output_index[1]) + std::to_string(output_index[2]) + "]"
  //    + "; Distance = " + std::to_string(output.first) + "\n";

  //  break;
  //}
  //default:
  //  break;
  //}

  //std::ofstream myFile(file_output);
  //myFile << output_message;
  //myFile.close();

  //std::cout << "Finished successfully\n";

  //return EXIT_SUCCESS;
}