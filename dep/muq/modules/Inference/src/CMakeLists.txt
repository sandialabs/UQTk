# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

CreateCompileGroup(
    INFERENCE_FILTERING
    "Inference components related to filtering, smoothing, and data assimilation."
    muqInference
    "MODELING_LINEARALGEBRA"
    "EIGEN3"
    ""
    Filtering/KalmanFilter.cpp
    Filtering/KalmanSmoother.cpp
)
