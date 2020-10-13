#ifndef ANYMAT_H_
#define ANYMAT_H_

#include <boost/any.hpp>

#include <typeinfo>
#include <typeindex>
#include <unordered_map>

namespace muq {
  namespace Modeling {


    class AnyVec{

    public:
        AnyVec(boost::any& objIn) : obj(objIn){};

        unsigned Size(const int dim=-1);
        double Norm();

        /* AnyMat operator*(AnyMat const& otherMat) const; */

        /* // Access operators */
        double  operator()(int i) const;
        double& operator()(int i);

        boost::any operator+(boost::any const& otherMat) const;
        boost::any operator+(AnyVec const& otherMat) const;

        boost::any operator-(boost::any const& otherMat) const;
        boost::any operator-(AnyVec const& otherMat) const;

        /* AnyMat& operator+=(AnyMat const& otherMat); */
        /* AnyMat& operator-=(AnyMat const& otherMat); */

        //boost::any& operator=(boost::any& otherMat);
        //boost::any& operator=(AnyMat otherMat);


        /* AnyMat Solve(AnyMat const& rhs); */

        /* AnyMat Cholesky(); */

        /* AnyMat Zero(int rows, int cols=-1) const; */
        /* AnyMat Ones(int rows, int cols=-1) const; */

        /* bool IsZero() const; */

        /* AnyMat Identity(int rows, int cols) const; */

        /* double LogDeterminant() const; */

        /* AnyMat Concatenate(AnyMat const& otherMat, int axis=-1) const; */


    private:
        boost::any& obj;

    };
  }
}
#endif
