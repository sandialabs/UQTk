#include <MUQ/Utilities/HDF5/H5Object.h>


int main()
{
    // Open the HDF5 file for reading and writing
    auto f = muq::Utilities::OpenFile("Data.h5");
        
    // Create a new group.  If it already exists, this will do nothing.
    f.CreateGroup("/NewGroup");
        
    // Add an attribute to the group.
    f["/NewGroup"].attrs["Some Metadata"] = "Created with MUQ!";
        
    // Store a vector in the group
    f["/NewGroup/Ones"] = Eigen::VectorXd::Ones(10);

    // Groups can also be stored and accessed.
    auto g = f["/NewGroup"];
    Eigen::VectorXd ones = g["/Ones"];
    std::cout << "\nThe content of /NewGroup/Ones is: \n" << ones.transpose() << std::endl;
    
    // Add some vector-valued metadata to the new dataset.
    f["/NewGroup/Ones"].attrs["Meta2"] = Eigen::VectorXd::Random(2).eval();
    
    // Store another vector in a different group.  The group is automatically generated
    f["/AnotherGroup/Zeros"] = Eigen::MatrixXd::Zero(5,5);

    // The content in a dataset can casted to an Eigen::MatrixX using the ".eval()" function.
    std::cout << "\nThe content of /AnotherGroup/Zeros is:\n" << f["/AnotherGroup/Zeros"].eval() << std::endl;

    return 0;
}
