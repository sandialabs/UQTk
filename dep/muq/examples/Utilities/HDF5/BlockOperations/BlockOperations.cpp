#include <MUQ/Utilities/HDF5/H5Object.h>


int main()
{
    // Open the HDF5 file for reading and writing
    auto f = muq::Utilities::OpenFile("Data.h5");
        
    // Store a vector in the group
    Eigen::VectorXd baseVector = Eigen::VectorXd::LinSpaced(11,0,10);
    f["/Vector"] = baseVector;

    // Now, extract the first two components of the dataset
    Eigen::VectorXd tempVector = f["/Vector"].head(2);

    std::cout << tempVector.transpose() << "\nvs\n"
	      << baseVector.head(2).transpose() << std::endl << std::endl;
    
    // Extract a length 3 segment from the middle of the dataset, staring at index 2.
    tempVector = f["/Vector"].segment(2,3);

    std::cout << tempVector.transpose() << "\nvs\n"
	      << baseVector.segment(2,3).transpose() << std::endl << std::endl;

    // Insert a vector into the middle of the dataset
    std::cout << "Old vector = "
	      << f["/Vector"].eval().transpose() << std::endl;
    
    f["/Vector"].segment(3,2) = Eigen::VectorXd::Zero(2).eval();
    
    std::cout << "New vector = "
	      << f["/Vector"].eval().transpose() << std::endl << std::endl;


    // Create a matrix dataset from the outer product of the [0,1,2...] vector
    f["/Matrix"] = (baseVector*baseVector.transpose()).eval();

    // Grab a block in the middle of matrix
    std::cout << "4x5 matrix block = \n" << f["/Matrix"].block(2,3,4,5).eval() << std::endl;

    // It is also possible to use all of the other Eigen block operations (e.g., .col(), .row(), .bottomLeftCorner(), .topRows(), .leftCols(), etc....)
    
    return 0;
}
