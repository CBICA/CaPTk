CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

SET( PROJECT_NAME Preprocessing )

CAPTK_ADD_PROJECT( ${PROJECT_NAME} 1.6.2 )

ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/ cbica_toolkit ) # takes care of find_itk

CAPTK_INITIAL_SETUP()

# update submodule
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/greedy/CMakeLists.txt" )
  EXECUTE_PROCESS(COMMAND git submodule update --init ${PROJECT_SOURCE_DIR}/greedy
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RESULT_VARIABLE git_result
  OUTPUT_VARIABLE git_ver)
ENDIF()

# ADD_SUBDIRECTORY(greedy)

### greedy is not picked up by macOS
SET(GREEDY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/greedy )
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${GREEDY_SOURCE_DIR}
  ${GREEDY_SOURCE_DIR}/src
  ${GREEDY_SOURCE_DIR}/src/ITKFilters/include
)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/greedy/src/GreedyVersion.cxx.in
  ${CMAKE_CURRENT_BINARY_DIR}/GreedyVersion.cxx @ONLY IMMEDIATE)
SET(GREEDY_HEADERS
  greedy/src/ITKFilters/include/FastLinearInterpolator.h
  greedy/src/ITKFilters/include/FastWarpCompositeImageFilter.h
  greedy/src/ITKFilters/include/FastWarpCompositeImageFilter.txx
  greedy/src/ITKFilters/include/JacobianDeterminantImageFilter.h
  greedy/src/ITKFilters/include/JacobianDeterminantImageFilter.txx
  greedy/src/ITKFilters/include/MultiComponentImageMetricBase.h
  greedy/src/ITKFilters/include/MultiComponentImageMetricBase.txx
  greedy/src/ITKFilters/include/MultiComponentMetricReport.h
  greedy/src/ITKFilters/include/MultiComponentMutualInfoImageMetric.h
  greedy/src/ITKFilters/include/MultiComponentMutualInfoImageMetric.txx
  greedy/src/ITKFilters/include/MultiComponentNCCImageMetric.h
  greedy/src/ITKFilters/include/MultiComponentNCCImageMetric.txx
  greedy/src/ITKFilters/include/MultiImageAffineMSDMetricFilter.h
  greedy/src/ITKFilters/include/MultiImageAffineMSDMetricFilter.txx
  greedy/src/ITKFilters/include/MultiImageOpticalFlowImageFilter.h
  greedy/src/ITKFilters/include/MultiImageOpticalFlowImageFilter.txx
  greedy/src/ITKFilters/include/OneDimensionalInPlaceAccumulateFilter.h
  greedy/src/ITKFilters/include/OneDimensionalInPlaceAccumulateFilter.txx
  greedy/src/ITKFilters/include/SimpleWarpImageFilter.h
  greedy/src/ITKFilters/include/SimpleWarpImageFilter.txx
  greedy/src/ITKFilters/include/itkGaussianInterpolateImageFunction.h
  greedy/src/ITKFilters/include/itkOptVectorLinearInterpolateImageFunction.h
  greedy/src/ITKFilters/include/itkOptVectorLinearInterpolateImageFunction.txx
  greedy/src/lddmm_common.h
  greedy/src/lddmm_data.h
  greedy/src/AffineCostFunctions.h
  greedy/src/AffineTransformUtilities.h
  greedy/src/GreedyAPI.h
  greedy/src/GreedyException.h
  greedy/src/GreedyParameters.h
  greedy/src/MultiImageRegistrationHelper.h
  greedy/src/CommandLineHelper.h
)

# Define greedy library files
SET(GREEDY_SRCS
  greedy/src/lddmm_data.cxx
  greedy/src/GreedyAPI.cxx
  greedy/src/GreedyParameters.cxx
  greedy/src/MultiImageRegistrationHelper.cxx
  greedy/src/AffineCostFunctions.cxx
  greedy/src/greedy_main.cxx
  ${CMAKE_CURRENT_BINARY_DIR}/GreedyVersion.cxx
)
CAPTK_ADD_EXECUTABLE( greedy "${GREEDY_HEADERS};${GREEDY_SRCS}" "${LIBNAME_CBICATK};${OpenMP_LIBRARIES}" )
###

INCLUDE_DIRECTORIES( ${CACHED_INCLUDE_DIRS} )

CAPTK_ADD_EXECUTABLE( ${PROJECT_NAME} "${APPLICATION_SOURCES};${APPLICATION_DEPENDS}" ${LIBNAME_CBICATK} )

ADD_DEPENDENCIES( ${PROJECT_NAME} ${LIBNAME_CBICATK} greedy )