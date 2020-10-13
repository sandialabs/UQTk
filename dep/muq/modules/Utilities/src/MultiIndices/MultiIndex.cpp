#include "MUQ/Utilities/MultiIndices/MultiIndex.h"

#include <iostream>
#include <stdexcept>

using namespace muq::Utilities;


MultiIndex::MultiIndex(unsigned lengthIn) : length(lengthIn),
                                            maxValue(0),
                                            totalOrder(0)
{}

MultiIndex::MultiIndex(unsigned lengthIn, unsigned val) : MultiIndex(lengthIn)
{
  for(int i=0; i<length; ++i){
    SetValue(i, val);
  }
}

MultiIndex::MultiIndex(Eigen::RowVectorXi const& indIn) : MultiIndex(indIn.size())
{
  maxValue = 0;
  totalOrder = 0;

  for(int i=0; i<indIn.size(); ++i){
    if( indIn[i] > 0 ){
      nzInds[i] = indIn[i];
      maxValue = std::max<int>(maxValue, indIn[i]);
      totalOrder += indIn[i];
    }
  }
}

MultiIndex::MultiIndex(std::initializer_list<unsigned> const& indIn) : MultiIndex(indIn.size())
{
  maxValue = 0;
  totalOrder = 0;

  unsigned i = 0;
  for(auto it = indIn.begin(); it != indIn.end(); ++it){
    if( *it > 0 ){
      nzInds[i] = *it;

      maxValue = std::max<int>(maxValue, *it);
      totalOrder += *it;

      i++;
    }
  }
}


Eigen::RowVectorXi MultiIndex::GetVector() const
{
  Eigen::RowVectorXi output = Eigen::RowVectorXi::Zero(length);

  for(auto it = nzInds.begin(); it!=nzInds.end(); ++it)
    output(it->first) = it->second;

  return output;
}

bool MultiIndex::SetValue(unsigned ind, unsigned val)
{
  if(ind>length){
    throw std::out_of_range("Tried to set the value of index " + std::to_string(ind) + " on an multiindex with only " + std::to_string(length) + " components.");
  }else{

    bool foundIndex;
    if (val > 0) {
      auto it = nzInds.find(ind);
      foundIndex = it!=nzInds.end();
      if(it != nzInds.end()){
        it->second = val;
      }else{
        nzInds[ind] = val;
      }
    } else {
      foundIndex = nzInds.erase(ind) > 0;
    }


    // Update the total and maximum order values after updating multi index
    totalOrder = 0;
    maxValue = 0;

    for (const auto& value : nzInds){
      totalOrder += value.second;
      maxValue = std::max(maxValue, value.second);
    }

    return foundIndex;
  }
}

unsigned int MultiIndex::NumNz() const
{
    unsigned int numNz = 0;
    for(auto& part : nzInds)
      numNz += int(part.second > 0);

    return numNz;
}


unsigned MultiIndex::MultiIndex::GetValue(unsigned ind) const
{
  if(ind>length){
    throw std::out_of_range("Tried to access index " + std::to_string(ind) + " of a multiindex with only " + std::to_string(length) + " components.");
  }else{
    auto searchIter = nzInds.find(ind);
    if(searchIter != nzInds.end()){
      return searchIter->second;
    }else{
      return 0;
    }
  }
}


void MultiIndex::SetLength(unsigned newLength)
{
  if(newLength > length){
    length = newLength;
  }else{

    auto it = nzInds.begin();
    while(it!= nzInds.end()){
      if (it->first >= newLength) {
	      it = nzInds.erase(it);
      } else {
	      it++;
      }
    }

    // Update the stored summaries
    length = newLength;
    maxValue = 0;
    totalOrder = 0;
    for(auto it = nzInds.begin(); it!=nzInds.end(); ++it){
      maxValue = std::max(maxValue, it->second);
      totalOrder += it->second;
    }

  }
}

bool MultiIndex::operator!=(const MultiIndex &b) const{

  if( (b.length != length) || (b.maxValue != maxValue) || (b.totalOrder != totalOrder))
    return true;

  if(b.nzInds.size() != nzInds.size())
    return true;

  // Loop through the nonzero indices
  auto bit = b.nzInds.begin();
  auto it = nzInds.begin();
  for(; it!=nzInds.end(); ++it){
    if(it->first != bit->first)
      return true;
    if(it->second != bit->second)
      return true;
  }

  return false;
}

bool MultiIndex::operator==(const MultiIndex &b) const{
  return !( *this != b);
}

bool MultiIndex::operator>(const MultiIndex &b) const{
  return b<(*this);
}

bool MultiIndex::operator<(const MultiIndex &b) const{

  if(totalOrder<b.totalOrder){
    return true;
  }else if(totalOrder>b.totalOrder){
    return false;
  }else if(maxValue<b.maxValue){
    return true;
  }else if(maxValue>b.maxValue){
    return false;
  }else{

    for(int i=0; i<std::min<unsigned>(length, b.length); ++i){
      if(GetValue(i)<b.GetValue(i)){
        return true;
      }else if(GetValue(i)>b.GetValue(i)){
        return false;
      }
    }

    // it should never get to this point unless the multiindices are equal
    return false;
  }

}

MultiIndex& MultiIndex::operator+=(const MultiIndex &b) {
  for(int i=0; i<length; ++i){
    SetValue(i, GetValue(i) + b.GetValue(i));
  }
  return *this;
}

MultiIndex& MultiIndex::operator++() {
  MultiIndex ones (this->GetLength(), 1);
  return (*this)+=ones;
}

MultiIndex MultiIndex::operator+(const MultiIndex &b) const{
  MultiIndex ret(*this);
  return ret += b;
}

MultiIndex& MultiIndex::operator-=(const MultiIndex &b) {
  for(int i=0; i<length; ++i){
    unsigned diff = 0;
    if (GetValue(i) > b.GetValue(i)) // Prevent "negative" unsigned result
      diff = GetValue(i) - b.GetValue(i);
    SetValue(i, diff);
  }
  return *this;
}

MultiIndex& MultiIndex::operator--() {
  MultiIndex ones (this->GetLength(), 1);
  return (*this)-=ones;
}

MultiIndex MultiIndex::operator-(const MultiIndex &b) const{
  MultiIndex ret(*this);
  return ret -= b;
}

std::string MultiIndex::ToString() const {
  std::string out;
  for(int i=0; i<GetLength(); ++i){
    if (i > 0)
      out += " ";
    out += std::to_string(GetValue(i));
  }
  return out;
}

std::ostream& muq::Utilities::operator<< (std::ostream &out, const MultiIndex &ind)
{
  out << ind.GetVector().transpose();
  return out;
}
