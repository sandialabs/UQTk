#include "MUQ/Utilities/HDF5/H5Object.h"
#include <iostream>

using namespace muq::Utilities;

void H5Object::DeepCopy(H5Object const& otherObj)
{
  file->Copy(path, otherObj.file, otherObj.path);
  children = otherObj.children;
  isDataset = otherObj.isDataset;
}

H5Object& H5Object::operator=(H5Object const& otherObj)
{
  DeepCopy(otherObj);
  return *this;
}

void H5Object::ExactCopy(H5Object const& otherObj)
{
  file = otherObj.file;
  attrs = otherObj.attrs;

  path = otherObj.path;
  children = otherObj.children;
  isDataset = otherObj.isDataset;
}

H5Object& H5Object::CreatePlaceholder(std::string const& grpName)
{
    auto pathParts = SplitString(grpName);

    if(pathParts.second.length()==0)
    {
      children[pathParts.first].ExactCopy( H5Object(file, path+grpName, true) );
      return children[pathParts.first];
    }

    if(children.find( pathParts.first ) == children.end()){
      if((path.length()==1) && (path.at(0)=='/'))
      {
	       file->CreateGroup(pathParts.first);
	       children[pathParts.first].ExactCopy( H5Object(file, pathParts.first ,false) );
      }
      else
      {
	       file->CreateGroup(path + pathParts.first);
	       children[pathParts.first].ExactCopy( H5Object(file, path + pathParts.first ,false) );
      }

      return children[pathParts.first].CreatePlaceholder(pathParts.second);
    }
    else
    {
      return children[pathParts.first].CreatePlaceholder(pathParts.second);
    }
};


H5Object& H5Object::CreateGroup(std::string const& grpName)
{
    // If the group already exists, just return it.
    if(file->DoesGroupExist(path + grpName))
	return (*this)[grpName];

    // Otherwise, recursively create the necessary groups.
    auto pathParts = SplitString(grpName);

    if(children.find( pathParts.first ) == children.end())
    {
	file->CreateGroup(path + pathParts.first);
	children[pathParts.first].ExactCopy( H5Object(file, path + pathParts.first, false) );

	if(pathParts.second.length()!=0)
	{
	    return children[pathParts.first].CreateGroup(pathParts.second);
	}
	else
	{
	    return children[pathParts.first];
	}
    }
    else
    {
	return children[pathParts.first].CreateGroup(pathParts.second);
    }
};




H5Object& H5Object::operator[](std::string const& targetPath)
{
  if(isDataset || (targetPath.length()==0))
  {
    return *this;
  }
  else
  {
    auto pathParts = SplitString(targetPath);

    if(children.find( pathParts.first ) != children.end())
    {
      return children.at(pathParts.first)[pathParts.second];
    }
    else
    {
      return CreatePlaceholder(targetPath);
    }
  }
};


// const H5Object& H5Object::operator[](std::string const& path) const
//   {
//       // If the full path is /this/is/the/path, then the split path will be
//       // /this and /is/the/path

//       if(isDataset || (path.length()==0))
//       {
// 	  return *this;
//       }
//       else
//       {
// 	  auto pathParts = SplitString(path);

// 	  return children.at(pathParts.first)[pathParts.second];
//       }
//   };

BlockDataset H5Object::block(unsigned startRow, unsigned startCol, unsigned numRows, unsigned numCols) const
{
    assert(isDataset);

    Eigen::VectorXi shape = file->GetDataSetSize(path);
    assert(shape.size()>0);
    assert(startRow+numRows <= shape(0));
    if(shape.size()>1){
      assert(startCol+numCols <= shape(1));
    }else{
      assert(numCols==1);
    }

    return BlockDataset(path, file, startRow, startCol, numRows, numCols);
}

BlockDataset H5Object::topLeftCorner(unsigned numRows, unsigned numCols) const
{
    return block(0,0,numRows,numCols);
}

BlockDataset H5Object::bottomLeftCorner(unsigned numRows, unsigned numCols) const
{
    return block(rows()-numRows,0,numRows,numCols);
}

BlockDataset H5Object::topRightCorner(unsigned numRows, unsigned numCols) const
{
    return block(0, cols()-numCols, numRows, numCols);
}

BlockDataset H5Object::bottomRightCorner(unsigned numRows, unsigned numCols) const
{
    return block(rows()-numRows, cols()-numCols, numRows, numCols);
}

BlockDataset H5Object::topRows(unsigned numRows) const
{
    return block(0, 0, numRows, cols());
}

BlockDataset H5Object::bottomRows(unsigned numRows) const
{
    return block(0, cols()-numRows, numRows, cols());
}

BlockDataset H5Object::leftCols(unsigned numCols) const
{
    return block(0, 0, rows(), numCols);
}

BlockDataset H5Object::rightCols(unsigned numCols) const
{
    return block(0,cols()-numCols, rows(), numCols);
}

BlockDataset H5Object::col(unsigned col) const
{
    assert(isDataset);
    Eigen::VectorXi shape = file->GetDataSetSize(path);


    return block(0,col,shape(0),1);
}

BlockDataset H5Object::row(unsigned row) const
{
    assert(isDataset);
    Eigen::VectorXi shape = file->GetDataSetSize(path);

    return block(row,0,1,shape(1));
}

BlockDataset H5Object::segment(unsigned startInd, unsigned numInds) const
{
    return block(startInd, 0, numInds, 1);
}

BlockDataset H5Object::head(unsigned numInds) const
{
    // Make sure the object is one dimensional
    Eigen::VectorXi shape = file->GetDataSetSize(path);

    if(shape.size()==0){
      std::cerr << "\nERROR: The dataset, " << path << ", does not exist.\n" << std::endl;
      assert(shape.size()>0);
    }

    if(shape.size()==1){
      return block(0,0,numInds,1);
    }else if(shape(0)==1){
      return block(0,0,1,numInds);
    }else if(shape(1)==1){
      return block(0,0,numInds,1);
    }else{
      std::cerr << "\nERROR: The head() function requires the dataset to be one dimensional and \"" << path << "\" does not seem to be one dimensional.\n" << std::endl;
      assert(false);
      return block(0,0,numInds,1);
    }
}

BlockDataset H5Object::tail(unsigned numInds) const
{
    return block(rows()-numInds, 0, numInds, 1);
}


unsigned H5Object::rows() const
{
    if(!isDataset)
      throw std::runtime_error("Attempted to call H5Object::rows() on a group object.");

    Eigen::VectorXi shape = file->GetDataSetSize(path);
    return shape(0);
}

unsigned H5Object::cols() const
{
  if(!isDataset)
    throw std::runtime_error("Attempted to call H5Object::cols() on a group object.");

    Eigen::VectorXi shape = file->GetDataSetSize(path);
    if(shape.size()==1)
      return 1;
    else
      return shape(1);
}

unsigned H5Object::size() const
{
  if(!isDataset)
    throw std::runtime_error("Attempted to call H5Object::size() on a group object.");

    Eigen::VectorXi shape = file->GetDataSetSize(path);
    return shape.prod();
}

double H5Object::operator()(int i) const
{
    if(isDataset)
    {
	return file->ReadPartialMatrix(path,i,0,1,1)(0);
    }
    else
    {
	assert(false);
    }
};
  
double H5Object::operator()(int i, int j) const
{
    if(isDataset)
    {
     	return file->ReadPartialMatrix(path, i,j,1,1)(0,0);
    }
    else
    {
	assert(false);
    }
};



void H5Object::Flush()
{
    file->FlushFile();
}
  

void H5Object::Print(std::string prefix) const
{
    std::cout << prefix + path << std::endl;
    for(auto& child : children)
	child.second.Print(prefix + "  ");
}


/// Recursively add children to create an HDF5 file hierarchy
H5Object muq::Utilities::AddChildren(std::shared_ptr<HDF5File>        file,
                                     std::string               const& groupName)
{

    if(file->IsDataSet(groupName))
	    return H5Object(file,groupName,true);

    // Set up the current object
    H5Object output(file, groupName, false);

    // Add the children
    std::vector<std::string> children = file->GetChildren(groupName);
    
    for(auto& childPath : children)
    {
	std::string fullChildPath = groupName;
	if(groupName.at(fullChildPath.length()-1)!='/')
	    fullChildPath += "/";
	fullChildPath += childPath;
	
	output.children[fullChildPath].ExactCopy( AddChildren(file, fullChildPath) );
    }
    return output;
}

/// Open an HDF5 file and return the root object
/**

 */
H5Object muq::Utilities::OpenFile(std::string const& filename)
{
    
    std::shared_ptr<HDF5File> file = std::make_shared<HDF5File>(filename);
    return AddChildren(file, "/");
}
