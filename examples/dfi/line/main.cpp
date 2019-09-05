#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <math.h>
#include <iostream>
//for setting stdout precision
#include <iomanip>

//#include "mcmc.h"
//#include "quad.h"
//#include "arrayio.h"
//#include "dsfmt_add.h"
//#include <sys/time.h>
//#include "arraytools.h"
#include <ctime>
#include <string>
#include <vector>
//#include "math.h"

#include "dfi.h"
//#include "molfrac.h"
//#include "def.h"

//#include "foo.h"

/*
#ifndef FOO_BAR_BAZ_H_
#define FOO_BAR_BAZ_H_
#endif  // FOO_BAR_BAZ_H_
*/

//seed
//uint32_t seed=13;

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
		//DFI dfiObject2("dfi.input");
        
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
        }/*
        //refit the consistent data 
        else if(strcmp(argv[1],"evidence")==0){
                cout<<"Computing model evidence"<<endl;
		DFI dfiObject1;
                dfiObject1.computeEvidence();
        }
	*/
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

//	DFI dfiObject1;

//	dfiObject1.dataInference();

//	DFI dfiObject2("foo");
//	exit(1);
	
	std::cout<<std::endl;
	std::cout<<"======================"<<std::endl;
	std::cout <<"EXITING DFI"<<std::endl;
	std::cout<<"======================"<<std::endl;
	
	return 0;
}

