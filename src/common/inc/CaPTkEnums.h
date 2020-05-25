#pragma once

namespace CAPTK
{
  enum ImageExtension
  {
    NIfTI = 0, DICOM
  };

  enum ImageModalityType
  {
    IMAGE_TYPE_UNDEFINED = 0, IMAGE_TYPE_T1, IMAGE_TYPE_T1CE, IMAGE_TYPE_T2,
    IMAGE_TYPE_T2FLAIR, IMAGE_TYPE_AX, IMAGE_TYPE_FA, IMAGE_TYPE_RAD, IMAGE_TYPE_TR,
    IMAGE_TYPE_PERFUSION, IMAGE_TYPE_DTI, IMAGE_TYPE_RECURRENCE_OUTPUT, IMAGE_TYPE_PP, IMAGE_TYPE_CT,
    IMAGE_TYPE_PET, IMAGE_TYPE_PSR, IMAGE_TYPE_PH, IMAGE_TYPE_RCBV, IMAGE_TYPE_SEG,
    IMAGE_TYPE_ATLAS, IMAGE_TYPE_PARAMS, IMAGE_TYPE_SUDOID, IMAGE_TYPE_NEAR, IMAGE_TYPE_FAR,
    IMAGE_MAMMOGRAM, IMAGE_TYPE_FEATURES
  };

  enum ClassificationConfigurationType
  {
    CONF_TYPE_UNDEFINED = 0, CONF_TYPE_KFOLD_CV, CONF_TYPE_DOUBLE, CONF_TYPE_SPLIT_TRAIN , CONF_TYPE_SPLIT_TEST
  };

  //! The modality strings that are used in the GUI 
  static const char ImageModalityString[CAPTK::ImageModalityType::IMAGE_TYPE_FEATURES + 1][15] =
  { "DEF", "T1", "T1Gd", "T2",
  "FLAIR", "DTI_AX", "DTI_FA", "DTI_RAD", "DTI_TR",
  "PERFUSION", "DTI", "REC", "PP", "CT",
  "PET", "pSR", "PH", "RCBV", "SEG",
  "ATLAS", "PARAMS", "SUDOID", "NEAR", "FAR",
  "FFDM", "FEAT" };

  enum MachineLearningApplicationSubtype
  {
    TRAINING = 0, TESTING
  };

  enum ApplicationCallingSVM
  {
    Default = 0, Recurrence, Survival, SBRT
  };

  enum SURVIVAL_SVM_PARAMS
  {
    BESTC_6MONTH = 32, BESTG_6MONTH = 1 / 32, BESTC_18MONTH = 1 / 2, BESTG_18MONTH = 1 / 32
  };

  enum GLISTR_OUTPUT_LABELS
  {
    EDEMA = 2, TUMOR = 4, VENT = 7, NONENHANCING = 1, ALL = 0
  };

  enum ORIGINAL_ATLAS_LABELS
  {
    VENT1 = 3, VENT2 = 8, VENT3 = 232, VENT4 = 233
  };

  enum MOLECULAR_SUBTYPES
  {
    PRONEURAL = 1, NEURAL = 2, MESSENCHYMAL = 3, CLASSICAL = 4
  };

  enum VOXEL_STATUS
  {
    OFF = 0, ON
  };

}