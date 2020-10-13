#ifndef FENICSPIECE_H
#define FENICSPIECE_H

#include "MUQ/Modeling/Dolfin/SwigExtract.h"
#include "MUQ/Modeling/WorkPiece.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include <dolfin/la/GenericVector.h>

#include <dolfin/fem/LinearVariationalProblem.h>
#include <dolfin/fem/LinearVariationalSolver.h>
#include <dolfin/fem/Form.h>

#include <dolfin/function/Function.h>
#include <dolfin/la/EigenVector.h>

#include <string>

#include <iostream>
#include <fstream>

#include <cxxabi.h>

std::string demangle_typename(const char* name) {

    int status = -4; // some arbitrary value to eliminate the compiler warning

    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(name, NULL, NULL, &status),
        std::free
    };

    return (status==0) ? res.get() : name ;
}


namespace muq{
    namespace Modeling{
        
        class FenicsPiece : public WorkPiece
        {
            
        public:
            FenicsPiece(pybind11::object              const& problemIn,
                        pybind11::object              const& outputFieldIn,
                        std::vector<pybind11::object> const& inputs) : WorkPiece(inputs.size(), 1){
                
                ExtractInputs(inputs);
                
                // First, make sure the input object is actually an instance of LinearVariationalProblem
                CheckProblemType(problemIn, "LinearVariationalProblem");
                CheckProblemType(outputFieldIn, "Function");
                
                problem = SwigExtract(problemIn);
                outputField = SwigExtract(outputFieldIn);
            };
            
            
            void EvaluateFunc(std::vector<pybind11::object> const& inputs)
            {
                // Set all of the expressions
                for(auto& expr : inputExprs)
                {
                    CheckProblemType(inputs.at(expr.first),"list");
                    
                    auto xList = pybind11::list(inputs.at(expr.first));
                    for(unsigned i=0; i<expr.second.second.size(); ++i)
                        expr.second.first.attr(expr.second.second.at(i).c_str())  = xList[i];
                }
                
                // Set all the functions
                for(auto& f : inputFuncs)
                {
                    CheckProblemType(inputs.at(f.first),"Function");
                    *f.second = *SwigExtract(inputs.at(f.first)).Cast<std::shared_ptr<dolfin::Function>>();
                }
                // Update each of the expressions
//        exprs.at(0).attr("bc_val") = leftVal;
//        exprs.at(1).attr("bc_val") = rightVal;
                
                //*inputField = *SwigExtract(f).Cast<std::shared_ptr<dolfin::Function>>();
                dolfin::LinearVariationalSolver solver(problem);
                solver.solve();
            }

            virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override
            {
                // Loop over the function inputs
                for(auto it = inputFuncs.begin(); it!= inputFuncs.end(); ++it)
                {
                    // This is a reference to the vector defining the Dolfin function that we want to update
                    Eigen::VectorXd& inVec = *std::dynamic_pointer_cast<dolfin::EigenVector>(it->second->vector())->vec();

                    // This is a reference to the vector containing the new values that we want to insert
                    Eigen::Ref<Eigen::Matrix<double, -1, -1, 1, -1, -1>, 0, Eigen::OuterStride<-1>> vec = boost::any_cast<Eigen::Ref<Eigen::Matrix<double, -1, -1, 1, -1, -1>, 0, Eigen::OuterStride<-1>> >(inputs.at(it->first).get());
                                            
                    // Set the new values and update the function
                    inVec = vec;
                    it->second->update();
                }


                // Loop over the expression inputs
                for(auto it = inputExprs.begin(); it!=inputExprs.end(); ++it)
                {
                    
                    unsigned inputInd = it->first;

                    std::pair<pybind11::object, std::vector<std::string>> &exprDef = it->second;
                    
                    // This is a reference to the vector containing the new values that we want to insert
                    if(inputs.at(inputInd).get().type() == typeid(double)){
                        double val = boost::any_cast<double>(inputs.at(inputInd).get());

                        assert(exprDef.second.size()==1);
                        
                        exprDef.first.attr(exprDef.second.at(0).c_str()) = pybind11::cast(val);
                        
                    }else{
                        Eigen::Ref<Eigen::Matrix<double, -1, -1, 1, -1, -1>, 0, Eigen::OuterStride<-1>> vec = boost::any_cast<Eigen::Ref<Eigen::Matrix<double, -1, -1, 1, -1, -1>, 0, Eigen::OuterStride<-1>> >(inputs.at(inputInd).get());
                        assert(exprDef.second.size()==vec.rows());
                        for(int k=0; k<exprDef.second.size(); ++k){
                            exprDef.first.attr(exprDef.second.at(k).c_str()) = pybind11::cast(vec(k,0));
                        }
                    }
                    
                }


                // Solve the system!
                dolfin::LinearVariationalSolver solver(problem);
                solver.solve();

                std::shared_ptr<dolfin::EigenVector> vec = std::dynamic_pointer_cast<dolfin::EigenVector>(outputField->vector());
                assert(vec);

                outputs.resize(1);
                outputs.at(0) = Eigen::VectorXd(*vec->vec());
            }
            
            boost::any EvaluateVec(Eigen::Ref<const Eigen::VectorXd> const& x, std::vector<double> vals)
            {
                //Eigen::VectorXd& inVec = std::dynamic_pointer_cast<dolfin::EigenVector>(inputField->vector())->vec();
                //assert(inVec);
                
                // Update the value of the input field
                //inVec = x;
                //inputField->update();
                
                // Update each of the expressions
                //exprs.at(0).attr("bc_val") = leftVal;
                //exprs.at(1).attr("bc_val") = rightVal;
                
                dolfin::LinearVariationalSolver solver(problem);
                solver.solve();
                
                std::shared_ptr<dolfin::EigenVector> vec = std::dynamic_pointer_cast<dolfin::EigenVector>(outputField->vector());
                assert(vec);
                
                return *vec->vec();
            };
            
        private:
            
            static void CheckProblemType(pybind11::object const& obj, std::string const& requiredType)
            {
                std::string typeName = pybind11::handle(obj).ptr()->ob_type->tp_name;
                
                if(typeName.compare(requiredType))
                {
                    throw std::invalid_argument("FenicsPiece constructor was given an instance of \"" + typeName + "\" but requires an instance of \"" + requiredType + "\"");
                }
            }
            
            void ExtractInputs(std::vector<pybind11::object> const& inputs){
                
                
                for(unsigned i=0; i<inputs.size(); ++ i){
                    
                    
                    // If the input is a list or tuple, then it should define an expression
                    if(pybind11::isinstance<pybind11::list>(inputs.at(i))){
                        
                        pybind11::list input(inputs[i]);
                        
                        pybind11::object expr = input[0];
                        CheckProblemType(expr, "CompiledExpression");
                        
                        CheckProblemType(input[1],"list");
                        pybind11::list part2(input[1]);
                        
                        if(pybind11::isinstance<pybind11::list>(part2) || pybind11::isinstance<pybind11::tuple>(part2)){
                            
                            std::vector<std::string> names;
                            for(auto& name : part2){
                                names.push_back(name.cast<std::string>());
                            }
                            inputExprs[i] = std::make_pair(expr, names);
                            
                        }
                        
                    }else{
                        CheckProblemType(inputs.at(i),"Function");
                        inputFuncs[i] = SwigExtract(inputs.at(i));
                    }
                }
            }
            
            std::map<unsigned, std::pair<pybind11::object, std::vector<std::string>>> inputExprs;
            std::map<unsigned, std::shared_ptr<dolfin::Function>> inputFuncs;
            
            std::shared_ptr<dolfin::Function> outputField;
            std::shared_ptr<dolfin::LinearVariationalProblem> problem;
        };
        
    }// namespace Modeling
}// namespace muq


#endif
