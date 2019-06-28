project (UQTk) 

# array tests
add_subdirectory (ArrayReadAndWrite) 
add_subdirectory (ArrayDelColumn   ) 
add_subdirectory (Array1DMiscTest  ) 
add_subdirectory (Array2DMiscTest  ) 
add_subdirectory (ArraySortTest    ) 

# # other tests
add_subdirectory (MultiIndexTest   ) 
add_subdirectory (CorrTest         )

# # quadrature class tests
add_subdirectory (QuadLUTest       ) 

# # MCMC class tests
add_subdirectory (MCMC2dTest       ) 
add_subdirectory (MCMCRandomTest   ) 
add_subdirectory (MCMCNestedTest   ) 

# # PCBasis and PCSet tests
add_subdirectory (Deriv1dTest      )
add_subdirectory (SecondDeriv1dTest) 
add_subdirectory (GradHessianTest  ) 
add_subdirectory (GradientPCETest  ) 
add_subdirectory (PCE1dTest        )
add_subdirectory (PCEImplTest      )
add_subdirectory (Hessian2dTest    ) 

# # BCS tests
add_subdirectory (BCS1dTest        ) 
add_subdirectory (BCS2dTest        ) 
