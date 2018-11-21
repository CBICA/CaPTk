/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#ifndef GREEDYAPI_H
#define GREEDYAPI_H

#include "GreedyParameters.h"
#include "GreedyException.h"
#include "lddmm_data.h"
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_random.h>
#include <map>

#include "itkCommand.h"

template <typename T, unsigned int V> class MultiImageOpticalFlowHelper;

namespace itk {
  template <typename T, unsigned int D1, unsigned int D2> class MatrixOffsetTransformBase;

}

/**
 * This is the top level class for the greedy software. It contains methods
 * for deformable and affine registration.
 */
template <unsigned int VDim, typename TReal = double>
class GreedyApproach
{
public:

  typedef GreedyApproach<VDim, TReal> Self;

  typedef LDDMMData<TReal, VDim> LDDMMType;
  typedef typename LDDMMType::ImageBaseType ImageBaseType;
  typedef typename LDDMMType::ImageType ImageType;
  typedef typename LDDMMType::ImagePointer ImagePointer;
  typedef typename LDDMMType::VectorImageType VectorImageType;
  typedef typename LDDMMType::VectorImagePointer VectorImagePointer;
  typedef typename LDDMMType::CompositeImageType CompositeImageType;
  typedef typename LDDMMType::CompositeImagePointer CompositeImagePointer;

  typedef vnl_vector_fixed<TReal, VDim> VecFx;
  typedef vnl_matrix_fixed<TReal, VDim, VDim> MatFx;

  typedef std::vector<std::vector<double> > MetricLogType;

  typedef MultiImageOpticalFlowHelper<TReal, VDim> OFHelperType;

  typedef itk::MatrixOffsetTransformBase<TReal, VDim, VDim> LinearTransformType;

  struct ImagePair {
    ImagePointer fixed, moving;
    VectorImagePointer grad_moving;
    double weight;
  };


  int Run(GreedyParameters &param);

  int RunDeformable(GreedyParameters &param);

  int RunAffine(GreedyParameters &param);

  int RunBrute(GreedyParameters &param);

  int RunReslice(GreedyParameters &param);

  int RunInvertWarp(GreedyParameters &param);

  int RunRootWarp(GreedyParameters &param);

  int RunAlignMoments(GreedyParameters &param);

  int RunJacobian(GreedyParameters &param);

  /**
   * Add an image that is already in memory to the internal cache, and
   * associate it with a filename. This provides a way for images already
   * loaded in memory to be passed in to the Greedy API while using the
   * standard parameter structures.
   *
   * Normally, images such as the fixed image are passed as part of the
   * GreedyParameters object as filenames. For example, we might set
   *
   *   param.inputs[0].fixed = "/tmp/goo.nii.gz";
   *
   * However, if we are linking to the greedy API from another program and
   * already have the fixed image in memory, we can use the cache mechanism
   * instead.
   *
   *   greedyapi.AddCachedInputObject("FIXED-0", myimage);
   *   param.inputs[0].fixed = "FIXED-0";
   *
   * The API will check the cache before loading the image. The type of the
   * object in the cache must match the type of the object expected internally,
   * which is VectorImage for most images. If not, an exception will be
   * thrown.
   *
   * Note that the cache does not use smart pointers to refer to the objects
   * so it's the caller's responsibility to keep the object pointed to while
   * the API is being used.
   */
  void AddCachedInputObject(std::string &string, itk::Object *object);

  /**
   * Get the metric log - values of metric per level. Can be called from
   * callback functions and observers
   */
  const MetricLogType &GetMetricLog() const;



protected:

  typedef std::map<std::string, itk::Object *> ImageCache;
  ImageCache m_ImageCache;

  // A log of metric values used during registration - so metric can be looked up
  // in the callbacks to RunAffine, etc.
  std::vector< std::vector<double> > m_MetricLog;

  // This function reads the image from disk, or from a memory location mapped to a
  // string. The first approach is used by the command-line interface, and the second
  // approach is used by the API, allowing images to be passed from other software
  template <class TImage>
  itk::SmartPointer<TImage> ReadImageViaCache(const std::string &filename);

  vnl_matrix<double> ReadAffineMatrixViaCache(const TransformSpec &ts);

  void WriteAffineMatrixViaCache(const std::string &filename, const vnl_matrix<double> &Qp);

  void ReadImages(GreedyParameters &param, OFHelperType &ofhelper);

  void ReadTransformChain(const std::vector<TransformSpec> &tran_chain,
                          ImageBaseType *ref_space,
                          VectorImagePointer &out_warp);

  static vnl_matrix<double> MapAffineToPhysicalRASSpace(
      OFHelperType &of_helper, int level,
      LinearTransformType *tran);

  static void MapPhysicalRASSpaceToAffine(
      OFHelperType &of_helper, int level,
      vnl_matrix<double> &Qp,
      LinearTransformType *tran);

  void RecordMetricValue(double val);

  // Compute the moments of a composite image (mean and covariance matrix of coordinate weighted by intensity)
  void ComputeImageMoments(CompositeImageType *image, const std::vector<double> &weights, VecFx &m1, MatFx &m2);

  class AbstractAffineCostFunction : public vnl_cost_function
  {
  public:

    AbstractAffineCostFunction(int n_unknowns) : vnl_cost_function(n_unknowns) {}
    virtual vnl_vector<double> GetCoefficients(LinearTransformType *tran) = 0;
    virtual void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran) = 0;
    virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g) = 0;
  };

  /**
   * Pure affine cost function - parameters are elements of N x N matrix M.
   * Transformation takes place in voxel coordinates - not physical coordinates (for speed)
   */
  class PureAffineCostFunction : public AbstractAffineCostFunction
  {
  public:

    typedef GreedyApproach<VDim, TReal> ParentType;

    // Construct the function
    PureAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);

    // Get the parameters for the specified initial transform
    vnl_vector<double> GetCoefficients(LinearTransformType *tran);

    // Get the transform for the specificed coefficients
    void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);

    // Get the preferred scaling for this function given image dimensions
    virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

    // Cost function computation
    virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

  protected:

    // Data needed to compute the cost function
    GreedyParameters *m_Param;
    OFHelperType *m_OFHelper;
    GreedyApproach<VDim, TReal> *m_Parent;
    bool m_Allocated;
    int m_Level;

    // Storage for the gradient of the similarity map
    VectorImagePointer m_Phi, m_GradMetric, m_GradMask;
    ImagePointer m_Metric, m_Mask;

    // Last set of coefficients evaluated
    vnl_vector<double> last_coeff;
  };

  /**
   * Physical space affine cost function - parameters are elements of affine transform in
   * physical RAS space.
   */
  class PhysicalSpaceAffineCostFunction : public AbstractAffineCostFunction
  {
  public:
    typedef GreedyApproach<VDim, TReal> ParentType;

    PhysicalSpaceAffineCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);
    virtual vnl_vector<double> GetCoefficients(LinearTransformType *tran);
    virtual void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);
    virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);
    virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

    void map_phys_to_vox(const vnl_vector<double> &x_phys, vnl_vector<double> &x_vox);

  protected:
    PureAffineCostFunction m_PureFunction;

    // Voxel to physical transforms for fixed, moving image
    typedef vnl_matrix_fixed<double, VDim, VDim> Mat;
    typedef vnl_vector_fixed<double, VDim> Vec;

    Mat Q_fix, Q_mov, Q_fix_inv, Q_mov_inv;
    Vec b_fix, b_mov, b_fix_inv, b_mov_inv;

    vnl_matrix<double> J_phys_vox;
  };

  /** Abstract scaling cost function - wraps around another cost function and provides scaling */
  class ScalingCostFunction : public AbstractAffineCostFunction
  {
  public:

    // Construct the function
    ScalingCostFunction(AbstractAffineCostFunction *pure_function, const vnl_vector<double> &scaling)
      : AbstractAffineCostFunction(pure_function->get_number_of_unknowns()),
        m_PureFunction(pure_function), m_Scaling(scaling) {}

    // Get the parameters for the specified initial transform
    vnl_vector<double> GetCoefficients(LinearTransformType *tran);

    // Get the transform for the specificed coefficients
    void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);

    // Cost function computation
    virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

    const vnl_vector<double> &GetScaling() { return m_Scaling; }

  protected:

    // Data needed to compute the cost function
    AbstractAffineCostFunction *m_PureFunction;
    vnl_vector<double> m_Scaling;
  };

  /** Cost function for rigid registration */
  class RigidCostFunction : public AbstractAffineCostFunction
  {
  public:
    typedef vnl_vector_fixed<double, VDim> Vec3;
    typedef vnl_matrix_fixed<double, VDim, VDim> Mat3;
    typedef GreedyApproach<VDim, TReal> ParentType;

    RigidCostFunction(GreedyParameters *param, ParentType *parent, int level, OFHelperType *helper);
    vnl_vector<double> GetCoefficients(LinearTransformType *tran);
    void GetTransform(const vnl_vector<double> &coeff, LinearTransformType *tran);
    virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g);

    // Get the preferred scaling for this function given image dimensions
    virtual vnl_vector<double> GetOptimalParameterScaling(const itk::Size<VDim> &image_dim);

    // Create a random set of parameters, such that on average point C_fixed maps to point C_mov
    vnl_vector<double> GetRandomCoeff(const vnl_vector<double> &xInit, vnl_random &randy, double sigma_angle, double sigma_xyz,
                                     const Vec3 &C_fixed, const Vec3 &C_moving);

  protected:

    Mat3 GetRotationMatrix(const Vec3 &q);
    Vec3 GetAxisAngle(const Mat3 &R);

    // We wrap around a physical space affine function, since rigid in physical space is not
    // the same as rigid in voxel space
    PhysicalSpaceAffineCostFunction m_AffineFn;

  };

  friend class GreedyApproach<VDim, TReal>::PureAffineCostFunction;

};

#endif // GREEDYAPI_H
