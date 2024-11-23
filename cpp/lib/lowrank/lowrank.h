/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file lowrank.h
/// \author Prashan Rai  2016 -
/// \brief Header file for Low Rank representation class

#ifndef CANLOWRANK_H_SEEN
#define CANLOWRANK_H_SEEN

#include "Array1D.h"
#include "Array2D.h"
#include "arraytools.h"



// Canonical tensor class definition
class CanonicalTensor {
private:
  Array1D<double> core_; // components are normalization constants of each rank one term
  Array1D< Array2D<double> > space_; // Each 2D array column contains univariate function coefficient. Number of arrays is equal to order of tensor
  Array1D<int> size_; // size of each mode
  int order_; // function dimension
public:
  /// \brief default constructor
  CanonicalTensor();
  /// \brief constructor
  CanonicalTensor(Array1D<double>& core, Array1D< Array2D<double> > space);
  /// \brief Destructor
  ~CanonicalTensor() {};
  /// \brief initialize as per str
  void create(string generator, int rank, Array1D<int>& size);
  /// \brief print the details of canonical tensor
  void display();
  /// \brief get the rank of canonical tensor
  size_t Rank(){return core_.Length();}
  /// \brief get order of the tensor
  size_t getOrder(){return order_;}
  /// \brief get the size of the tensor
  Array1D<int> getSize(){return size_;}
  /// \brief get the core
  Array1D<double> getCore(){return core_;}
  /// \brief get the core element (in this case, normalization coefficents of rank one terms)
  double getCoreVal(int rank){return core_(rank);}
  /// \brief set the core
  void setCore(Array1D<double>& core){core_ = core;}
  /// \brief set the core element at rank r to alpha
  void setCoreVal(int r, double alpha){core_(r) = alpha;}
  /// \brief Display space along a given dimension
 void displaySpace(int dim){printarray(space_[dim]);};
  /// \brief Get all univariate coefficients of a given dimension
  Array2D<double> getSpace(int dim){return space_[dim];};
  /// \brief Get rank univariate coefficients of a given dimension
  Array1D<double> getSpaceVec(int rank, int dim);
  /// \brief Set the univariate function coefficients of a given rank in given dim
  void setSpaceVec(int rank, int dim, Array1D<double>& column);
  /// \brief Print a canonical tensor
  void printTensor();
  /// \brief add two canonical tensors (to be done later)
  static CanonicalTensor add(CanonicalTensor& u, CanonicalTensor& v);
  /// \brief norm
  double norm();
  /// \brief Calculate scalar product of two canonical tensors
  static double prodscal(CanonicalTensor& u, CanonicalTensor& v);
};

class FunctionalBases{
public:
  /// \brief Default Constructor
  FunctionalBases() {};
  /// \brief Destructor
  ~FunctionalBases() {};
  /// \brief Get the number of dimensions for a functional bases type
  virtual size_t getNumOrder() = 0;
  /// \brief Evaluate bases
  virtual void functionEval(Array2D<double>& xval, int dim, Array2D<double>& bb) = 0;// dummy
};

class PCBases: public FunctionalBases{
private:
  Array1D<string> pctype_;
  Array1D<int> order_;
public:
  /// \brief Constructor
  PCBases(){};
  /// \brief Destructor
  ~PCBases(){};
  /// \brief Constructor that sets private data members
  PCBases(Array1D<string> pctype, Array1D<int> order);
  /// \brief Returns the number of dimensions with PC bases
  size_t getNumOrder(){return order_.Length();};
  /// \brief Get the order of pc bases in a particular dimension
  size_t getOrder(int dim){return order_(dim);};
  /// \brief Get the pc basis type in a particular dimension
  string getPCType(int dim){return pctype_(dim);}
  /// \brief Evaluate the basis functions in dimension dim for all data points and store in bb
  void functionEval(Array2D<double>& xval, int dim, Array2D<double>& bb);
};

class PLBases: public FunctionalBases{
private:
  Array1D<int> order_;
public:
  /// \brief Constructor
  PLBases(){};
  /// \brief Destructor
  ~PLBases(){};
  /// \brief Constructor that sets private data members
  PLBases(Array1D<int> order){order_ = order;};
  /// \brief Returns the number of dimensions with polynomial bases
  size_t getNumOrder(){return order_.Length();};
  /// \brief Get the order of a particular dimension
  size_t getOrder(int dim){return order_(dim);};
  /// \brief Evaluate the basis functions in dimension dim for all data points and store in bb
  void functionEval(Array2D<double>& xval, int dim, Array2D<double>& bb);
};

class CanonicalTensorALSLeastSquares{
private:
  /// \brief bases contains the type of bases
  Array1D<FunctionalBases* > bases_;
  /// \brief basesEval contains the evaluation of basis functions at sample points
  Array1D< Array2D<double> > basesEval_;
  /// \brief rank
  int rank_;
  /// \brief algorithm type
  string alg_;
  /// \brief tolerance to stop ALS
  double tol_;
  /// \brief initialization type
  string initializationType_;
  /// \brief max iterations in ALS
  double maxIterations_;
  /// \brief ragularization parameter
  double lambda_;
public:
  /// \brief default constructor
  CanonicalTensorALSLeastSquares();
  /// \brief Destructor
  ~CanonicalTensorALSLeastSquares() {};
  /// \brief Setting the bases as one of the derived class of functionalBases
  void setBases(Array1D<FunctionalBases* > bases){bases_ = bases;};
  /// \brief Setting basesEval as 2D arrays of basis functions evaluated at sample points
  void setBasesEval(Array1D <Array2D<double> >& basesEval){this->basesEval_ = basesEval;};
  /// \brief Set rank
  void setRank(int rank){rank_ = rank;};
    /// \brief Set number of iterations
  void setIteration(int maxiter){maxIterations_ = maxiter;};
  /// \brief Initialize using ones, zeros, random or greedy
  void setLambda(int lambda){lambda_ = lambda;};
  /// \brief Initialize using ones, zeros, random or greedy
  void initialization(CanonicalTensor& u);
  /// \brief ALS solver
  void solve(CanonicalTensor& u);
  /// \brief Greedy ALS Solver
  void solveGreedy(CanonicalTensor& u);
  /// \brief Direct Solver
  CanonicalTensor solveDirect(Array2D<double>& xval, Array1D<double>& yval,bool vflag);
  /// \brief Evaluate the basis function
  //template <typename T> Array1D <Array2D<double> > EvalBases(Array2D<double>& xx, T& BClass);
};

class FunctionalTensor{
private:
  CanonicalTensor tensor_;
  Array1D<FunctionalBases* > bases_;
  size_t fdims_;

public:
  /// \brief default constructor
  FunctionalTensor();
  /// \brief Constructors with inputs
  FunctionalTensor(CanonicalTensor& f);
  FunctionalTensor(CanonicalTensor& f, Array1D<FunctionalBases* >& bases);
  /// \brief Destructor
  ~FunctionalTensor(){};
  /// \brief Get the tensor
  CanonicalTensor getTensor(){return tensor_;};
  /// \brief Evaluate the tensor when we have FunctionalBases
  Array1D<double> tensorEval(Array2D<double>& xval);
  /// \brief Evaluate the tensor when we have bases evaluated at N samples
  Array1D<double> tensorEval(Array1D< Array2D<double> >& H, int N);
};

#endif /* CANLOWRANK_H_SEEN */


