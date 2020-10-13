#include "MUQ/Utilities/HDF5/Attributes.h"

muq::Utilities::Attribute& muq::Utilities::Attribute::operator=(std::string const& val)
{
    assert(file);
    file->WriteStringAttribute(path, name, val);
    return *this;
}



muq::Utilities::Attribute::operator std::string() const
{
    assert(file);
    return file->GetStringAttribute(path,name);
};


muq::Utilities::Attribute&  muq::Utilities::AttributeList::operator[](std::string const& attrName)
{
    assert(file);
    
    // If the attribute doesn't exist, create it and return the reference
	if( attributes.find(attrName) == attributes.end() )
	    attributes[attrName] = Attribute(file, path, attrName);
	
	return attributes[attrName];
};
