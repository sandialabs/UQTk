#include "MUQ/Utilities/HDF5/BlockDataset.h"

using namespace muq::Utilities;


BlockDataset& BlockDataset::operator=(boost::any const& val) {

    AnyWriterMapType& map = *GetAnyWriterMap();
    auto iter = map.find(val.type());
    if(iter == map.end()){
        std::cerr << "ERROR: MUQ does not know how to write a boost::any with underlying type \"" << val.type().name() << "\".  Currently implemented types are:\n";
        for(auto mapIter = map.begin(); mapIter!=map.end(); ++mapIter)
            std::cerr << "    " << mapIter->first.name() << std::endl;
        std::cerr << std::endl;

        assert(iter != map.end());
    }
    
    map[val.type()](val, *this);
    
    return *this;
}

std::shared_ptr<BlockDataset::AnyWriterMapType> BlockDataset::GetAnyWriterMap() {
    
    static std::shared_ptr<BlockDataset::AnyWriterMapType> map;

  if( !map )
      map = std::make_shared<BlockDataset::AnyWriterMapType>();

  return map;
}

REGISTER_HDF5BLOCK_ANYTYPE(double, double)
REGISTER_HDF5BLOCK_ANYTYPE(float, float)
REGISTER_HDF5BLOCK_ANYTYPE(int, int)
REGISTER_HDF5BLOCK_ANYTYPE(unsigned, unsigned)
REGISTER_HDF5BLOCK_ANYTYPE(MatrixXd, Eigen::MatrixXd)
REGISTER_HDF5BLOCK_ANYTYPE(MatrixXi, Eigen::MatrixXi)
REGISTER_HDF5BLOCK_ANYTYPE(MatrixXf, Eigen::MatrixXf)
REGISTER_HDF5BLOCK_ANYTYPE(VectorXf, Eigen::VectorXf)
REGISTER_HDF5BLOCK_ANYTYPE(VectorXd, Eigen::VectorXd)
REGISTER_HDF5BLOCK_ANYTYPE(VectorXi, Eigen::VectorXi)

