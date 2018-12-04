#ifndef _NiftiDataManager_h_
#define _NiftiDataManager_h_

//#include "CAPTk.h"
#include "cbicaITKSafeImageIO.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"


using ImageTypeFloat4D = itk::Image< float, 4 >;
using VectorVectorDouble = std::vector< std::vector < double > >;
using VectorDouble = std::vector < double >;
class NiftiDataManager
{
public:
  //!Constructor
  NiftiDataManager();
  
  //!Destructor
  ~NiftiDataManager();

  /**
  \brief Loads training data from given input images and a label image
  \param labelImagePointer Segmentated image
  \param recurrenceMaskImage Near region drawing mask
  \param nonrecurrenceMaskImage Far region drawing mask
  \param t1ceImagePointer Pointer of T1-Gd image
  \param t2flairImagePointer Pointer of T2 FLAIR image
  \param t1ImagePointer Pointer of T1 image
  \param t2ImagePointer Pointer of T2 image
  \param perfImagePointerNifti Pointer of Perfusion image in NIfTI format
  \param perfImagePointerDicom Pointer of Perfusion image in DICOM format
  \param axImagePointer   Pointer of Axial diffusivity image
  \param faImagePointer   Pointer of fractional anisotropy image
  \param radImagePointer  Pointer of radial diffusivity image
  \param trImagePointer   Pointer of trace image
  \param pNearIntensities Perfusion features of near points
  \param pFarIntensities  Perfusion features of far points
  \param mNearIntensities Other features (structural + diffusion) of near points
  \param mFarIntensities  Other features (structural + diffusion) of far points
  \param tNearIntensities Distance feature of near points
  \param tFarIntensities  Distance features of far points
  \param imagetype      Tells whether image is NIfTI or DICOM
  \param useT1Data      Tells whether to use T1 feature or not
  \param useT2Data      Tells whether to use T2 feature or not
  \param useT1CEData    Tells whether to use T1CE feature or not
  \param useT2FlairData Tells whether to use T2 Flair feature or not
  \param useDTIData     Tells whether to use DTI features or not
  \param usePerfData    Tells whether to use perfusion feature or not
  \param useDistData    Tells whether to use distance feature or not
  */
  void LoadTrainingData(ImageTypeFloat3D::Pointer labelImagePointer,
    ImageTypeFloat3D::Pointer recurrenceMaskImage,
    ImageTypeFloat3D::Pointer nonrecurrenceMaskImage,
    ImageTypeFloat3D::Pointer t1ceImagePointer,
    ImageTypeFloat3D::Pointer t2flairImagePointer,
    ImageTypeFloat3D::Pointer t1ImagePointer,
    ImageTypeFloat3D::Pointer t2ImagePointer,
    ImageTypeFloat4D::Pointer perfImagePointerNifti,
    ImageTypeFloat3D::Pointer axImagePointer,
    ImageTypeFloat3D::Pointer faImagePointer,
    ImageTypeFloat3D::Pointer radImagePointer,
    ImageTypeFloat3D::Pointer trImagePointer,
    VectorVectorDouble & pNearIntensities,
    VectorVectorDouble & pFarIntensities,
    VectorVectorDouble & mNearIntensities,
    VectorVectorDouble & mFarIntensities,
    VectorDouble & tNearIntensities,
    VectorDouble & tFarIntensities,
    int imagetype,
	bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData);


  template<class ImageType>
  typename ImageType::Pointer GetDistanceMap(const typename ImageType::Pointer &labelImagePointer);

  /**
  \brief Loads training data from given input images and a label image
  \param labelImagePointer Segmentated image
  \param nearpoints A vector having indices of near points
  \param farpoints  A vector having indices of far points
  \param t1ceImagePointer     Pointer of T1-Gd image
  \param t2flairImagePointer  Pointer of T2 FLAIR image
  \param t1ImagePointer       Pointer of T1 image
  \param t2ImagePointer       Pointer of T2 image
  \param perfImagePointer     Pointer of Perfusion image 
  \param dtiImagePointer      A vector of pointers of DTI images
  \param pNearIntensities Perfusion features of near points
  \param pFarIntensities  Perfusion features of far points
  \param mNearIntensities Other features (structural + diffusion) of near points
  \param mFarIntensities  Other features (structural + diffusion) of far points
  \param tNearIntensities Distance feature of near points
  \param tFarIntensities  Distance features of far points
  \param imagetype      Tells whether image is NIfTI or DICOM
  \param useT1Data      Tells whether to use T1 feature or not
  \param useT2Data      Tells whether to use T2 feature or not
  \param useT1CEData    Tells whether to use T1CE feature or not
  \param useT2FlairData Tells whether to use T2 Flair feature or not
  \param useDTIData     Tells whether to use DTI features or not
  \param usePerfData    Tells whether to use perfusion feature or not
  \param useDistData    Tells whether to use distance feature or not
  */
  void LoadTrainingData(ImageTypeFloat3D::Pointer labelImagePointer,
    VectorVectorDouble nearPoints,
    VectorVectorDouble farPoints,
    ImageTypeFloat3D::Pointer t1ceImagePointer,
    ImageTypeFloat3D::Pointer t2flairImagePointer,
    ImageTypeFloat3D::Pointer t1ImagePointer,
    ImageTypeFloat3D::Pointer t2ImagePointer,
    std::vector<ImageTypeFloat3D::Pointer> perfImagePointer,
    std::vector<ImageTypeFloat3D::Pointer> dtiImagePointer,
    VectorVectorDouble & pNearIntensities,
    VectorVectorDouble & pFarIntensities,
    VectorVectorDouble & mNearIntensities,
    VectorVectorDouble & mFarIntensities,
    VectorDouble & tNearIntensities,
    VectorDouble & tFarIntensities,
    int imagetype,
	bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData);

  /**
  \brief Loads test data from given test images and a label image

  \param t1ceImagePointer     Pointer of T1-Gd image
  \param t2flairImagePointer  Pointer of T2 FLAIR image
  \param t1ImagePointer       Pointer of T1 image
  \param t2ImagePointer       Pointer of T2 image
  \param perfImagePointerNifti Pointer of Perfusion image in NIfTI format
  \param perfImagePointerDicom Pointer of Perfusion image in DICOM format
  \param axImagePointer   Pointer of Axial diffusivity image
  \param faImagePointer   Pointer of fractional anisotropy image
  \param radImagePointer  Pointer of radial diffusivity image
  \param trImagePointer   Pointer of trace image
  \param labelImagePointer    Overall segmentated image
  \param dilatedEdemaPointer  Edema segmentated image
  \param pIntensities Perfusion features of test points
  \param mIntensities Other features (structural + diffusion) of test points
  \param tIntensities Distance feature of test points
  \param imagetype      Tells whether image is NIfTI or DICOM
  \param useT1Data      Tells whether to use T1 feature or not
  \param useT2Data      Tells whether to use T2 feature or not
  \param useT1CEData    Tells whether to use T1CE feature or not
  \param useT2FlairData Tells whether to use T2 Flair feature or not
  \param useDTIData     Tells whether to use DTI features or not
  \param usePerfData    Tells whether to use perfusion feature or not
  \param useDistData    Tells whether to use distance feature or not
  */
  std::vector<ImageTypeFloat3D::IndexType>  LoadTestData(ImageTypeFloat3D::Pointer t1ceImagePointer,
    ImageTypeFloat3D::Pointer t2flairImagePointer,
    ImageTypeFloat3D::Pointer t1ImagePointer,
    ImageTypeFloat3D::Pointer t2ImagePointer,
    ImageTypeFloat4D::Pointer perfImagePointerNifti,
    ImageTypeFloat3D::Pointer axImagePointer,
    ImageTypeFloat3D::Pointer faImagePointer,
    ImageTypeFloat3D::Pointer radImagePointer,
    ImageTypeFloat3D::Pointer trImagePointer,
    ImageTypeFloat3D::Pointer labelImagePointer,
    ImageTypeFloat3D::Pointer dilatedEdemaPointer,
    VectorVectorDouble & pIntensities,
    VectorVectorDouble & mIntensities,
    VectorDouble & tIntensities,
    int imagetype,
    bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData);

  /**
  \brief Loads test data from given test images and a label image
  \param t1ceImagePointer     Pointer of T1-Gd image
  \param t2flairImagePointer  Pointer of T2 FLAIR image
  \param t1ImagePointer       Pointer of T1 image
  \param t2ImagePointer       Pointer of T2 image
  \param perfImagePointer     Pointer of Perfusion image
  \param dtiImagePointer      A vector of pointers of DTI images
  \param labelImagePointer    Overall segmentated image
  \param dilatedEdemaPointer  Edema segmentated image
  \param pIntensities Perfusion features of test points
  \param mIntensities Other features (structural + diffusion) of test points
  \param tIntensities Distance feature of test points
  \param imagetype      Tells whether image is NIfTI or DICOM
  \param useT1Data      Tells whether to use T1 feature or not
  \param useT2Data      Tells whether to use T2 feature or not
  \param useT1CEData    Tells whether to use T1CE feature or not
  \param useT2FlairData Tells whether to use T2 Flair feature or not
  \param useDTIData     Tells whether to use DTI features or not
  \param usePerfData    Tells whether to use perfusion feature or not
  \param useDistData    Tells whether to use distance feature or not
  */

  std::vector<ImageTypeFloat3D::IndexType> LoadTestData(ImageTypeFloat3D::Pointer t1ceImagePointer,
    ImageTypeFloat3D::Pointer t2flairImagePointer,
    ImageTypeFloat3D::Pointer t1ImagePointer,
    ImageTypeFloat3D::Pointer t2ImagePointer,
    std::vector<ImageTypeFloat3D::Pointer> perfImagePointer,
    std::vector<ImageTypeFloat3D::Pointer> dtiImagePointer,
    ImageTypeFloat3D::Pointer labelImagePointer,
    ImageTypeFloat3D::Pointer dilatedEdemaPointer,
    VectorVectorDouble & pIntensities,
    VectorVectorDouble & mIntensities,
    VectorDouble & tIntensities,
    int imagetype,
    bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData);

  /**
  \brief Reads 3D Nifti image
  \param filename Input image file name
  */

  itk::Image<float, 3>::Pointer ReadNiftiImage(std::string filename);

  /**
  \brief Reads 4D Nifti image
  \param filename Input image file name
  */
  itk::Image<float, 4>::Pointer Read4DNiftiImage(std::string filename);

  /**
  \brief Reads image with specified dimensions and input pixel type
  \param filename Input image file name
  */
  template<class InputPixelType = float, class OutputPixelType = float, unsigned int VImageDimension = 3>
  typename itk::Image<OutputPixelType, VImageDimension>::Pointer  ReadImageWithDimAndInputPixelType(std::string filename);

};

template<class InputPixelType, class OutputPixelType, unsigned int VImageDimension>
typename itk::Image<OutputPixelType, VImageDimension>::Pointer  NiftiDataManager::ReadImageWithDimAndInputPixelType(std::string filename)
{
  return cbica::ReadImage< itk::Image< OutputPixelType, VImageDimension > >(filename);
}
template<class ImageType>
typename ImageType::Pointer NiftiDataManager::GetDistanceMap(const typename ImageType::Pointer &labelImagePointer)
{
	typedef itk::Image<unsigned char, 3>  UnsignedCharImageType;
	typedef itk::CastImageFilter< ImageType, UnsignedCharImageType > CastFilterType;
	typename CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(labelImagePointer);
	typedef  itk::SignedMaurerDistanceMapImageFilter< UnsignedCharImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;
	typename  SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();
	distanceMapImageFilter->SetInput(castFilter->GetOutput());
	distanceMapImageFilter->Update();
	typename ImageType::Pointer DistanceMap = distanceMapImageFilter->GetOutput();
	return DistanceMap;
}
#endif
