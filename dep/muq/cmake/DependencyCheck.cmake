
# The inference library requires the optimization library
if(Inference_build)
  set(Optimization_build ON)
endif(Inference_build)

# If incremental approximation is needed, make sure nlopt is available
IF(NOT MUQ_USE_NLOPT AND IncrementalApproximation_build)
	message(WARNING "IncrementalApproximation explicitly requires NLOPT.  NLOPT was not found, so IncrementalApproximation will not be compiled.")
	set(IncrementalApproximation_build OFF)
	set(IncrementalApproximation_tests OFF)
endif()


# If there is no gtest, turn off all tests
if(NOT MUQ_USE_GTEST)
  set(Inference_tests OFF)
  set(Optimization_tests OFF)
  set(Approximation_tests OFF)
  set(Geostats_tests OFF)
  set(UtilitiesAndModelling_tests OFF)
  set(PolynomialChaos_tests OFF)
  message(WARNING "MUQ_USE_GTEST is OFF. Thus, all tests has been turned off.")
endif(NOT MUQ_USE_GTEST)
