/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <math.h>
#include <iostream>
//for setting stdout precision
#include <iomanip>

#include <ctime>
#include <string>
#include <vector>

#include "dfi.h"

int main (int argc, char *argv[]) {

	std::cout<<"======================"<<std::endl;
	std::cout <<"WELCOME TO DFI"<<std::endl;
	std::cout<<"======================"<<std::endl;
	std::cout<<std::endl;

	//============================================
	if (argc==1){
                std::cout<<"Error: no input keyword supplied"<<std::endl;
	}
	//run a data inference
        else if (strcmp(argv[1],"datainference")==0) {
                std::cout<<"Running data inference"<<std::endl;
		DFI dfiObject1;
		//DFI dfiObject1("dfi.input");

	        dfiObject1.buildSurrogateModel();
	        dfiObject1.loadSurrogateModel();
	        dfiObject1.dataInference();
        }
        //refit the consistent data
        else if(strcmp(argv[1],"datarefit")==0){
                std::cout<<"Performing data refitting"<<std::endl;
		DFI dfiObject1;
	        dfiObject1.buildSurrogateModel();
	        //dfiObject1.loadSurrogateModel();
                dfiObject1.dataRefit();
        }
	//build surrogate model
        else if(strcmp(argv[1],"buildSurrogate")==0){
                cout<<"Building surrogate model"<<endl;
		DFI dfiObject1;
                dfiObject1.buildSurrogateModel();
	}
	//load surrogate model
        else if(strcmp(argv[1],"loadSurrogate")==0){
                cout<<"Loading surrogate model"<<endl;
		DFI dfiObject1;
                dfiObject1.loadSurrogateModel();
	}
	//test surrogate model against detailed model
        else if(strcmp(argv[1],"testSurrogate")==0){
                cout<<"Testing surrogate model"<<endl;
		DFI dfiObject1;
		dfiObject1.loadSurrogateModel();
		dfiObject1.testSurrogateModel();
	}
	//build KDE posteriors and pool
        else if(strcmp(argv[1],"buildKDE")==0){

		//check if necessary addtional arguments were supplied
		if (argc==2){
                	std::cout<<"Error: no dimensions selected for KDE"<<std::endl;
		}else{
                	cout<<"Building KDE posteriors"<<endl;

			//parse command line arguments to determine dimensions to use for building KDE density
			stringstream KDEparse;
			int parsedint;
			//std::cout<<KDEparse.str()<<std::endl;
			Array1D<int> KDEdimensions;
			for (int i=2;i<argc;i++){
				KDEparse.clear();
				KDEparse.str("");
				KDEparse.str(argv[i]);
				KDEparse>>parsedint;
				KDEdimensions.PushBack(parsedint);
			}
			std::cout<<"Will construct KDE posterior using MCMC chain dimension(s): ";
			for (int i=0;i<KDEdimensions.XSize()-1;i++){
				std::cout<<KDEdimensions(i)<<", ";
			}
			std::cout<<KDEdimensions(KDEdimensions.XSize()-1)<<std::endl;

			//check that all parsed integers are positive, >0
			bool parseCheck=true;
			for (int i=0;i<KDEdimensions.XSize();i++){
				if(KDEdimensions(i)<=0){
				parseCheck=false;
				}
			}

			//if all is well, proceed
			if (parseCheck){
				DFI dfiObject1;
                		dfiObject1.buildKDE(KDEdimensions);
			}else{
                		std::cout<<"Error: an input dimension is <=0. Exiting..."<<std::endl;
			}
        	}
	}
	//argument not found
        else{
                std::cout<<"Error: input keyword not recognized"<<std::endl;
        }
	//============================================

	std::cout<<std::endl;
	std::cout<<"======================"<<std::endl;
	std::cout <<"EXITING DFI"<<std::endl;
	std::cout<<"======================"<<std::endl;

	return 0;
}
