/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: SimpleWarpImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-10-29 11:19:00 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __MultiImageSimpleWarpImageFilter_h
#define __MultiImageSimpleWarpImageFilter_h
#include "itkImageBase.h"
#include "itkImageFunction.h"
#include "itkImageToImageFilter.h"
#include "itkPoint.h"
#include "itkFixedArray.h"
#include "itkVectorImage.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkInPlaceImageFilter.h"

/**
 * This class is used to perform mean square intensity difference type
 * registration with multiple images. The filter is designed for speed
 * of interpolation.
 */
template <class TFloat, unsigned int VDim>
class MultiImageOpticalFlowHelper
{
public:

  typedef itk::VectorImage<TFloat, VDim> MultiComponentImageType;
  typedef itk::Image<TFloat, VDim> FloatImageType;
  typedef itk::CovariantVector<TFloat, VDim> VectorType;
  typedef itk::Image<VectorType, VDim> VectorImageType;
  typedef itk::ImageBase<VDim> ImageBaseType;

  typedef std::vector<int> PyramidFactorsType;

  typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> LinearTransformType;

  /** Set default (power of two) pyramid factors */
  void SetDefaultPyramidFactors(int n_levels);

  /** Set the pyramid factors - for multi-resolution (e.g., 8,4,2) */
  void SetPyramidFactors(const PyramidFactorsType &factors);

  /** Add a pair of multi-component images to the class - same weight for each component */
  void AddImagePair(MultiComponentImageType *fixed, MultiComponentImageType *moving, double weight);

  /** Compute the composite image - must be run before any sampling is done */
  void BuildCompositeImages();

  /** Get the reference image for level k */
  ImageBaseType *GetReferenceSpace(int level);

  /** Get the reference image for level k */
  ImageBaseType *GetMovingReferenceSpace(int level);

  /** Perform interpolation - compute [(I - J(Tx)) GradJ(Tx)] */
  double ComputeOpticalFlowField(int level, VectorImageType *def, VectorImageType *result,
                                 double result_scaling = 1.0);

  double ComputeAffineMatchAndGradient(int level, LinearTransformType *tran,
                                       LinearTransformType *grad = NULL);


protected:

  // Pyramid factors
  PyramidFactorsType m_PyramidFactors;

  // Weights
  std::vector<double> m_Weights;

  // Vector of images
  typedef std::vector<typename MultiComponentImageType::Pointer> MultiCompImageSet;

  // Fixed and moving images
  MultiCompImageSet m_Fixed, m_Moving;

  // Composite image at each resolution level
  MultiCompImageSet m_FixedComposite, m_MovingComposite;

  void PlaceIntoComposite(FloatImageType *src, MultiComponentImageType *target, int offset);
  void PlaceIntoComposite(VectorImageType *src, MultiComponentImageType *target, int offset);
};

namespace itk
{

template<class TFloat, class TFloatArr, unsigned int VDim>
static void flatten_affine_transform(
    const MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform,
    TFloatArr *flat_array)
{
  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    flat_array[pos++] = transform->GetOffset()[i];
    for(int j = 0; j < VDim; j++)
      flat_array[pos++] = transform->GetMatrix()(i,j);
    }
}

template<class TFloat, class TFloatArr, unsigned int VDim>
static void unflatten_affine_transform(
   const TFloatArr *flat_array,
   MatrixOffsetTransformBase<TFloat, VDim, VDim> *transform,
   double scaling = 1.0)
{
  typename MatrixOffsetTransformBase<TFloat, VDim, VDim>::MatrixType matrix;
  typename MatrixOffsetTransformBase<TFloat, VDim, VDim>::OffsetType offset;

  int pos = 0;
  for(int i = 0; i < VDim; i++)
    {
    offset[i] = flat_array[pos++] * scaling;
    for(int j = 0; j < VDim; j++)
      matrix(i, j) = flat_array[pos++] * scaling;
    }

  transform->SetMatrix(matrix);
  transform->SetOffset(offset);
}


template<class TInputImage, class TOutputImage, class TDeformationImage>
class MultiImageOpticalFlowWarpTraits
{
public:
  typedef TDeformationImage TransformType;

  typedef typename TInputImage::InternalPixelType InputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputPixelType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TDeformationImage::ImageDimension );

  static DataObject *AsDataObject(TransformType *t) { return t; }
  static ImageBase<ImageDimension> *AsImageBase(TransformType *t) { return t; }

  static int GetResultAccumSize(int) { return 1; }

  static int GetStride(int) { return 1; }

  static bool InterpolateOutsideOverlapRegion() { return true; }

  static void TransformIndex(const itk::Index<ImageDimension> &pos,
                            TransformType *transform, long offset,
                            float *ptran)
    {
    typename TDeformationImage::InternalPixelType &def = transform->GetBufferPointer()[offset];
    for(int i = 0; i < ImageDimension; i++)
      ptran[i] = pos[i] + def[i];
    }

  static void PostInterpolate(
      const itk::Index<ImageDimension> &pos,
      const InputPixelType *pFix, const InputPixelType *pMov, int nComp,
      float *weight, float mask, double *summary, OutputPixelType &vOut)
  {
    for(int i = 0; i < ImageDimension; i++)
      vOut[i] = 0;

    const InputPixelType *pMovEnd = pMov + nComp;
    while(pMov < pMovEnd)
      {
      double del = (*pFix++) - *(pMov++);
      double delw = (*weight++) * del;
      for(int i = 0; i < ImageDimension; i++)
        vOut[i] += delw * *(pMov++);
      *summary += delw * del;
      }
  }
};

template<class TInputImage, class TOutputImage>
class MultiImageOpticalFlowAffineGradientTraits
{
public:

  typedef typename TInputImage::InternalPixelType InputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputPixelType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  typedef MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> TransformType;


  static DataObject *AsDataObject(TransformType *t) { return NULL; }
  static ImageBase<ImageDimension> *AsImageBase(TransformType *t) { return NULL; }

  static int GetResultAccumSize(int nComp) { return 1 + ImageDimension * (1 + ImageDimension); }

  static int GetStride(int) { return 1; }

  static bool InterpolateOutsideOverlapRegion() { return true; }

  static void TransformIndex(const itk::Index<ImageDimension> &pos,
                            TransformType *transform, long offset,
                            float *ptran)
    {
    for(int i = 0; i < ImageDimension; i++)
      {
      ptran[i] = transform->GetOffset()[i];
      for(int j = 0; j < ImageDimension; j++)
        ptran[i] += transform->GetMatrix()(i,j) * pos[j];
      }
    }

  static void PostInterpolate(
      const itk::Index<ImageDimension> &pos,
      const InputPixelType *pFix, const InputPixelType *pMov, int nComp,
      float *weight, float mask, double *summary, OutputPixelType &vOut)
  {
    const InputPixelType *pMovEnd = pMov + nComp;
    for(int i = 0; i < ImageDimension; i++)
      vOut[i] = 0;

    if(mask == 1.0)
      {
      while(pMov < pMovEnd)
        {
        double del = (*pFix++) - *(pMov++);
        double delw = (*weight++) * del;
        for(int i = 0; i < ImageDimension; i++)
          vOut[i] += delw * *(pMov++);
        *summary += delw * del;
        }

      for(int i = 0; i < ImageDimension; i++)
        {
        *(++summary) += vOut[i];
        for(int j = 0; j < ImageDimension; j++)
          {
          *(++summary) += vOut[i] * pos[j];
          }
        }
      }
    else if(mask > 0.0)
      {
      while(pMov < pMovEnd)
        {
        double del = (*pFix++) - *(pMov++);
        double delw = (*weight++) * del;
        //for(int i = 0; i < ImageDimension; i++)
        //  vOut[i] += delw * ( *(pMov++) * mask + del *
        //*summary += delw * del;
        }

      for(int i = 0; i < ImageDimension; i++)
        {
        *(++summary) += vOut[i];
        for(int j = 0; j < ImageDimension; j++)
          {
          *(++summary) += vOut[i] * pos[j];
          }
        }
      }


    /*
      */


  }
};




template<class TInputImage, class TOutputImage>
class MultiImageOpticalFlowAffineObjectiveTraits
{
public:
  typedef MultiImageOpticalFlowAffineGradientTraits<TInputImage,TOutputImage> SourceTraits;
  typedef typename SourceTraits::TransformType TransformType;

  typedef typename TInputImage::InternalPixelType InputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputPixelType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  static DataObject *AsDataObject(TransformType *t) { return NULL; }
  static ImageBase<ImageDimension> *AsImageBase(TransformType *t) { return NULL; }

  // We keep track of the average difference between fixed and interpolating moving
  // images as well as the size of the overlap region (i.e., number of voxels where
  // the measurement was obtained
  static int GetResultAccumSize(int) { return 2; }

  static int GetStride(int) { return 1 + ImageDimension; }

  static bool InterpolateOutsideOverlapRegion() { return true; }

  static void TransformIndex(const itk::Index<ImageDimension> &pos,
                            TransformType *transform, long offset,
                            float *ptran)
    {
    return SourceTraits::TransformIndex(pos, transform, offset, ptran);
    }

  static void PostInterpolate(
      const itk::Index<ImageDimension> &pos,
      const InputPixelType *pFix, const InputPixelType *pMov, int nComp,
      float *weight, float mask, double *summary, OutputPixelType &vOut)
  {
<<<<<<< HEAD
    double wdiff = 0.0;

    if(mask > 0.0)
      {
      for(int i = 0; i < nComp; i+=(1+ImageDimension))
        {
        double del = (*pFix++) - *(pMov++);
        double delw = (*weight++) * del;
        wdiff += delw * del;
        }

      summary[0] += wdiff * mask;
      summary[1] += mask;
=======
    double avgsqdiff;
    for(int i = 0; i < nComp; i+=(1+ImageDimension))
      {
      double del = (*pFix++) - *(pMov++);
      double delw = (*weight++) * del;
      avgsqdiff += delw * del;
>>>>>>> 43ef7496f075d647d8e516d0c8c81fc86f04a1ae
      }
    summary[0] += avgsqdiff;
    summary[1] += 1.0;
  }
};



/** \class MultiImageOpticalFlowImageFilter
 * \brief Warps an image using an input deformation field (for LDDMM)
 *
 * This filter efficiently computes the optical flow field between a
 * set of image pairs, given a transformation phi. This filter is the
 * workhorse of deformable and affine rigid registration algorithms that
 * use the mean squared difference metric. Given a set of fixed images F_i
 * and moving images M_i, it computes
 *
 *   v(x) = Sum_i w_i \[ F_i(x) - M_i(Phi(x)) ] \Grad M_i (Phi(x))
 *
 * The efficiency of this filter comes from combining the interpolation of
 * all the M and GradM terms in one loop, so that all possible computations
 * are reused
 *
 * The fixed and moving images must be passed in to the filter in the form
 * of VectorImages of size K and (VDim+K), respectively - i.e., the moving
 * images and their gradients are packed together.
 *
 * The output should be an image of CovariantVector type
 *
 * \warning This filter assumes that the input type, output type
 * and deformation field type all have the same number of dimensions.
 *
 */
template <class TInputImage, class TOutputImage, class TTransformTraits>
class ITK_EXPORT MultiImageOpticalFlowImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef MultiImageOpticalFlowImageFilter             Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiImageOpticalFlowImageFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension );

  /** Typedef to describe the output image region type. */
  typedef typename TInputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef TInputImage                                 InputImageType;
  typedef typename TInputImage::PixelType             InputPixelType;
  typedef typename TInputImage::InternalPixelType     InputComponentType;
  typedef TOutputImage                                OutputImageType;
  typedef typename OutputImageType::PixelType         OutputPixelType;
  typedef typename OutputPixelType::ComponentType     OutputComponentType;
  typedef typename OutputImageType::IndexType         IndexType;
  typedef typename OutputImageType::IndexValueType    IndexValueType;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::SpacingType       SpacingType;
  typedef typename OutputImageType::DirectionType     DirectionType;

  /** Information from the parent class */
  typedef typename TTransformTraits::TransformType    TransformType;
  typedef typename TransformType::Pointer             TransformPointer;

  /** Weight vector */
  typedef vnl_vector<float>                           WeightVectorType;
  typedef vnl_vector<double>                          SummaryType;


  /** typedef for base image type at the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

  /** Set the fixed image(s) */
  void SetFixedImage(InputImageType *fixed)
    { this->ProcessObject::SetInput("Primary", fixed); }

  /** Set the moving image(s) and their gradients */
  void SetMovingImageAndGradient(InputImageType *moving)
    { this->ProcessObject::SetInput("moving", moving); }

  /** Set the weight vector */
  itkSetMacro(Weights, WeightVectorType)
  itkGetConstMacro(Weights, WeightVectorType)

  /** Set the transform field. */
  void SetTransform(TransformType *transform)
    {
    m_Transform = transform;
    if(TTransformTraits::AsDataObject(transform))
      this->ProcessObject::SetInput("transform", TTransformTraits::AsDataObject(transform));
    }

  /** Summary results after running the filter */
  itkGetConstMacro(SummaryResult, SummaryType)

  /** This filter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information according the OutputSpacing, OutputOrigin
   * and the deformation field's LargestPossibleRegion. */
  virtual void GenerateOutputInformation();

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after
   * multi-threading. */
  virtual void AfterThreadedGenerateData();

protected:
  MultiImageOpticalFlowImageFilter();
  ~MultiImageOpticalFlowImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** SimpleWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            ThreadIdType threadId );

  void VerifyInputInformation() {}

  // Object to assist specializaiton
  struct DispatchBase {};
  template <unsigned int VDim> struct Dispatch : public DispatchBase {};

  /** Fast interpolation method */
  /*
  double OpticalFlowFastInterpolate(const Dispatch<3> &dispatch,
                                    float *cix,
                                    const InputComponentType *fixed_ptr,
                                    const InputComponentType *moving_ptr,
                                    OutputPixelType &outVector,
                                    int *movSize,
                                    int nComp,
                                    const InputComponentType *def_value);

  // Dummy implementation
  double OpticalFlowFastInterpolate(const DispatchBase &base,
                                    float *cix,
                                    const InputComponentType *fixed_ptr,
                                    const InputComponentType *moving_ptr,
                                    OutputPixelType &outVector,
                                    int *movSize,
                                    int nComp,
                                    const InputComponentType *def_value)
    { return 0.0; }
    */

  bool OpticalFlowFastInterpolate(const Dispatch<3> &dispatch,
                                  const InputComponentType *moving_ptr,
                                  int nComp, int stride, int *movSize,
                                  const InputComponentType *def_value,
                                  float *cix,
                                  InputComponentType *out);

  bool OpticalFlowFastInterpolate(const DispatchBase &base,
                                  const InputComponentType *moving_ptr,
                                  int nComp, int stride, int *movSize,
                                  const InputComponentType *def_value,
                                  float *cix,
                                  InputComponentType *out) { return true; }

  void OpticalFlowFastInterpolateWithMask(const Dispatch<3> &dispatch,
                                  const InputComponentType *moving_ptr,
                                  int nComp, int stride, int *movSize,
                                  const InputComponentType *def_value,
                                  float *cix,
                                  InputComponentType *out, float &outMask);

  void OpticalFlowFastInterpolateWithMask(const DispatchBase &base,
                                  const InputComponentType *moving_ptr,
                                  int nComp, int stride, int *movSize,
                                  const InputComponentType *def_value,
                                  float *cix,
                                  InputComponentType *out, float &outMask) {  }

private:
  MultiImageOpticalFlowImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Weight vector
  WeightVectorType                m_Weights;

  // Transform pointer
  TransformPointer                m_Transform;

  // Vector of accumulated data (difference, gradient of affine transform, etc)
  SummaryType                     m_SummaryResult;
  std::vector<SummaryType>        m_SummaryResultPerThread;
};






/**
 * This filter computes the similarity between a set of moving images and a
 * set of fixed images in a highly optimized way
 */
template <class TInputImage>
class ITK_EXPORT MultiImageAffineMSDMetricFilter :
    public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  /** Standard class typedefs. */
  typedef MultiImageAffineMSDMetricFilter              Self;
  typedef InPlaceImageFilter<TInputImage>              Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiImageAffineMSDMetricFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  /** Typedef to describe the output image region type. */
  typedef typename TInputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef TInputImage                                 InputImageType;
  typedef ImageBase<ImageDimension>                   ImageBaseType;
  typedef typename TInputImage::PixelType             InputPixelType;
  typedef typename TInputImage::InternalPixelType     InputComponentType;
  typedef typename InputImageType::IndexType          IndexType;
  typedef typename InputImageType::IndexValueType     IndexValueType;
  typedef typename InputImageType::SizeType           SizeType;
  typedef typename InputImageType::SpacingType        SpacingType;
  typedef typename InputImageType::DirectionType      DirectionType;

  /** Information from the parent class */
  typedef MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> TransformType;
  typedef typename TransformType::Pointer             TransformPointer;

  /** Weight vector */
  typedef vnl_vector<float>                           WeightVectorType;

  /** Set the fixed image(s) */
  void SetFixedImage(InputImageType *fixed)
    { this->ProcessObject::SetInput("Primary", fixed); }

  /** Set the moving image(s) and their gradients */
  void SetMovingImageAndGradient(InputImageType *moving)
    { this->ProcessObject::SetInput("moving", moving); }

  /** Set the weight vector */
  itkSetMacro(Weights, WeightVectorType)
  itkGetConstMacro(Weights, WeightVectorType)

  /** Whether to compute gradient */
  itkSetMacro(ComputeGradient, bool)
  itkGetConstMacro(ComputeGradient, bool)

  /** Set the transform field. */
  void SetTransform(TransformType *transform)
    { m_Transform = transform; }

  itkGetConstMacro(Transform, TransformType *)

  /** Value of the similarity objective after running the filter */
  itkGetConstMacro(MetricValue, double)

  /** The gradient (in the form of a transform) after running the filter */
  itkGetConstMacro(MetricGradient, TransformType *)



protected:
  MultiImageAffineMSDMetricFilter() : m_ComputeGradient(false) {}
  ~MultiImageAffineMSDMetricFilter() {}

  void PrintSelf(std::ostream& os, Indent indent) const
    { this->PrintSelf(os, indent); }

  /** SimpleWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for
   * ThreadedGenerateData(). */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            ThreadIdType threadId );

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void GenerateInputRequestedRegion();

  /** Override since input passed to output */
  virtual void EnlargeOutputRequestedRegion(DataObject *data);

  /** This method is used to set the state of the filter before
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after
   * multi-threading. */
  virtual void AfterThreadedGenerateData();

  /** Allocate outputs - just pass through the input */
  virtual void AllocateOutputs();

  void VerifyInputInformation() {}

  // Object to assist specializaiton
  struct DispatchBase {};
  template <unsigned int VDim> struct Dispatch : public DispatchBase {};

private:
  MultiImageAffineMSDMetricFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Weight vector
  WeightVectorType                m_Weights;

  // Transform pointer
  TransformPointer                m_Transform;

  // Whether the gradient is computed
  bool                            m_ComputeGradient;

  // Data accumulated for each thread
  struct ThreadData {
    double metric, mask;
    vnl_vector<double> gradient, grad_mask;
    ThreadData() : metric(0.0), mask(0.0),
      gradient(ImageDimension * (ImageDimension+1), 0.0),
      grad_mask(ImageDimension * (ImageDimension+1), 0.0) {}
  };

  std::vector<ThreadData>         m_ThreadData;

  // Vector of accumulated data (difference, gradient of affine transform, etc)
  double                          m_MetricValue;

  // Gradient
  TransformPointer                m_MetricGradient;
};









/** \class MultiImageOpticalFlowImageFilter
 * \brief Warps an image using an input deformation field (for LDDMM)
 *
 * This filter efficiently computes the optical flow field between a
 * set of image pairs, given a transformation phi. This filter is the
 * workhorse of deformable and affine rigid registration algorithms that
 * use the mean squared difference metric. Given a set of fixed images F_i
 * and moving images M_i, it computes
 *
 *   v(x) = Sum_i w_i \[ F_i(x) - M_i(Phi(x)) ] \Grad M_i (Phi(x))
 *
 * The efficiency of this filter comes from combining the interpolation of
 * all the M and GradM terms in one loop, so that all possible computations
 * are reused
 *
 * The fixed and moving images must be passed in to the filter in the form
 * of VectorImages of size K and (VDim+K), respectively - i.e., the moving
 * images and their gradients are packed together.
 *
 * The output should be an image of CovariantVector type
 *
 * \warning This filter assumes that the input type, output type
 * and deformation field type all have the same number of dimensions.
 *
 */
#ifdef SHAHAHA
template <class TInputImage, class TOutputImage, class TDeformationField = TOutputImage>
class ITK_EXPORT MultiImageOpticalFlowImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef MultiImageOpticalFlowImageFilter             Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiImageOpticalFlowImageFilter, ImageToImageFilter )

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension );

  /** Typedef to describe the output image region type. */
  typedef typename TInputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef TInputImage                                 InputImageType;
  typedef typename TInputImage::PixelType             InputPixelType;
  typedef typename TInputImage::InternalPixelType     InputComponentType;
  typedef TOutputImage                                OutputImageType;
  typedef typename OutputImageType::PixelType         OutputPixelType;
  typedef typename OutputPixelType::ComponentType     OutputComponentType;
  typedef typename OutputImageType::IndexType         IndexType;
  typedef typename OutputImageType::IndexValueType    IndexValueType;
  typedef typename OutputImageType::SizeType          SizeType;
  typedef typename OutputImageType::SpacingType       SpacingType;
  typedef typename OutputImageType::DirectionType     DirectionType;

  typedef itk::MatrixOffsetTransformBase<double, ImageDimension, ImageDimension> TransformType;

  /** Weight vector */
  typedef vnl_vector<float>                           WeightVectorType;

  /** typedef for base image type at the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

  /** Deformation field typedef support. */
  typedef TDeformationField                        DeformationFieldType;
  typedef typename DeformationFieldType::Pointer   DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType DisplacementType;

  /** Set the fixed image(s) */
  void SetFixedImage(InputImageType *fixed)
    { this->ProcessObject::SetInput("Primary", fixed); }

  /** Set the moving image(s) and their gradients */
  void SetMovingImageAndGradient(InputImageType *moving)
    { this->ProcessObject::SetInput("moving", moving); }

  /** Set the weight vector */
  itkSetMacro(Weights, WeightVectorType)
  itkGetConstMacro(Weights, WeightVectorType)

  /** Set the deformation field. */
  void SetDeformationField(DeformationFieldType *field)
    { this->ProcessObject::SetInput("deformation", field); }

  /** Set the affine transform - currently mutually exclusive with the deformation */
  void SetLinearTransform(TransformType *transform);

  /** Set constant scaling factor for the deformation field */
  itkSetMacro(DeformationScaling, float)
  itkGetConstMacro(DeformationScaling, float)

  /** Get the total energy of optical flow - only after Update has been called */
  itkGetConstMacro(TotalEnergy, double)

  /** This filter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information according the OutputSpacing, OutputOrigin
   * and the deformation field's LargestPossibleRegion. */
  virtual void GenerateOutputInformation();

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before 
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after 
   * multi-threading. */
  virtual void AfterThreadedGenerateData();

protected:
  MultiImageOpticalFlowImageFilter();
  ~MultiImageOpticalFlowImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** SimpleWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for 
   * ThreadedGenerateData(). */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            ThreadIdType threadId );

  void VerifyInputInformation() {}

  // Object to assist specializaiton
  struct DispatchBase {};
  template <unsigned int VDim> struct Dispatch : public DispatchBase {};

  /** Fast interpolation method */
  double OpticalFlowFastInterpolate(const Dispatch<3> &dispatch,
                                    float *cix,
                                    const InputComponentType *fixed_ptr,
                                    const InputComponentType *moving_ptr,
                                    OutputPixelType &outVector,
                                    int *movSize,
                                    int nComp,
                                    const InputComponentType *def_value);

  // Dummy implementation
  double OpticalFlowFastInterpolate(const DispatchBase &base,
                                    float *cix,
                                    const InputComponentType *fixed_ptr,
                                    const InputComponentType *moving_ptr,
                                    OutputPixelType &outVector,
                                    int *movSize,
                                    int nComp,
                                    const InputComponentType *def_value)
    { return 0.0; }

private:
  MultiImageOpticalFlowImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Scaling for the deformation field
  float                           m_DeformationScaling;

  // Weight vector
  WeightVectorType                m_Weights;

  // Linear transform
  typename TransformType::Pointer m_Transform;

  // Total energy - Sum |I_k - J_k|^2
  double                          m_TotalEnergy;
  std::vector<double>             m_TotalEnergyPerThread;

  // Gradient of the affine transform - computed when the deformation is null
  typename TransformType::Pointer m_GradTransform;
  std::vector<typename TransformType::Pointer> m_GradTransformPerThread;
};
#endif

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "MultiImageSimpleWarpImageFilter.txx"
#endif

#endif
