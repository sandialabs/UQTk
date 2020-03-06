/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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


#define MAX(a,b)        (((a) > (b)) ? (a) : (b))
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))


RefPtr<XMLElement> readXMLTree(const string inFileName);

void readXMLChainInput(RefPtr<XMLElement> xmlTree,MCMC* pmchain, Array1D<double>& chstart, int* nsteps, Array1D<int>& chainParamInd, Array1D<string>& priortype, Array1D<double>& priorparam1, Array1D<double>& priorparam2);
void readXMLModelInput(RefPtr<XMLElement> xmlTree,Array1D<double>& modelparams,Array1D<string>& modelparamnames,Array1D<double>& modelauxparams);
void readXMLDataInput(RefPtr<XMLElement> xmlTree, Array2D<double>& data,Array1D<double>& postparams, string* noiseType );
void readXMLUncInput(RefPtr<XMLElement> xmlTree, Array2D<double>& allPCcoefs,Array1D<int>& uncParamInd,int* outOrder, string* PCtype);


RefPtr<XMLElement> readXMLTree(const string inFileName)
{

// Compute the name of the input and write-back XML file
  unsigned long int xmlLPos = inFileName.find(".xml");
  unsigned long int xmlUPos = inFileName.find(".XML");

  string outFileExt = ".parsed.xml";

  string outFileName = inFileName;
  if(xmlLPos != string::npos){
    outFileName.replace(xmlLPos,4,outFileExt);
  }else if(xmlUPos != string::npos){
    outFileName.replace(xmlUPos,4,outFileExt);
  }else{
    throw Tantrum(std::string("The XML input file ") + inFileName + 
                  " does not have a valid (.xml or .XML) extension");
  }

  ifstream xmlInFile(inFileName.c_str(),ios::in);
  if(!xmlInFile){
    throw Tantrum(std::string("Could not open XML input file ") + inFileName);
  }
  ofstream xmlOutFile(outFileName.c_str(),ios::out);
  if(!xmlOutFile){
    throw Tantrum(std::string("Could not open XML output file ") + outFileName);
  }

  // Parse the input xml file
  RefPtr<XMLParser> parser = new XMLExpatParser;
   RefPtr<XMLElement> xmlTree = parser->parse(xmlInFile);

  // Write the parsed tree back to the write-back file
  dump_xml_tree(xmlTree, "", xmlOutFile);

return xmlTree;

}



void readXMLDataInput(RefPtr<XMLElement> xmlTree, Array2D<double>& data,Array1D<double>& postparams, string* noiseType )
{


  ///////////////////////////////////////////////////

  // Select the tree element with the input parameters for this run
  RefPtr<XMLElement> dataInput = xmlTree->get_child("Data");
  string filename = dataInput->attributes()->get("file","no_name_specified");
  *noiseType = dataInput->attributes()->get("noise","const_std");

  read_datafileVS(data,filename.c_str());
  RefPtr<XMLAttributeList> postpar = dataInput->get_child("noise_types")->get_child(*noiseType)->attributes();

  int npostparams=postpar->size();
  postparams.Resize(npostparams);
  postparams(0)=postpar->get_double("stdparam",0.1);

  return;
}


void readXMLModelInput(RefPtr<XMLElement> xmlTree, Array1D<double>& modelparams,Array1D<string>& modelparamnames,Array1D<double>& modelauxparams )
{


  RefPtr<XMLElement> params = xmlTree->get_child("ModelParameters");
  int np=params->count_children();
  for (int i=0;i<np;i++){
    modelparams.PushBack(params->get_child(i)->attributes()->get_double("value",0.e0));
    modelparamnames.PushBack(params->get_child(i)->label());
  }
  
// modelauxparams is a placeholder and currently not used


  return;
}

void readXMLUncInput(RefPtr<XMLElement> xmlTree, Array2D<double>& allPCcoefs,Array1D<int>& uncParamInd,int* outOrder, string* PCtype)
{

  RefPtr<XMLElement> pcInput = xmlTree->get_child("PC");
  *PCtype = pcInput->attributes()->get("type","GH");
  *outOrder = pcInput->attributes()->get_int("order",1);
  string inpc_type= pcInput->attributes()->get("inpc","default");

  RefPtr<XMLElement> inpc_types = pcInput->get_child("inpc_types");
  RefPtr<XMLElement> pt = inpc_types->get_child(inpc_type);  

  bool noFileFlag=true;
  string inpcfile;
  if(!strcmp(inpc_type.c_str(),"file")) {
    noFileFlag=false;
    inpcfile=pt->attributes()->get("name","input_pc.dat");

  }
  

 RefPtr<XMLElement> params = xmlTree->get_child("ModelParameters");
  int np=params->count_children();
  
  uncParamInd.Clear();
  Array1D<string> types(np);
  // count the number of uncertain params
  int ic=0;
  
  for (int i=0;i<np;i++){
    types(i)=params->get_child(i)->attributes()->get("type","no_type_given");
    if(!strcmp(types(i).c_str(),"uncertain")){
      uncParamInd.PushBack(i);
      ic++;
    }
  }
  int uncdim=ic;
  int npc=computeNPCTerms(uncdim,(*outOrder));
  allPCcoefs.Resize(npc,np,0.e0);
  for (int i=0;i<np;i++){
    
	double val=params->get_child(i)->attributes()->get_double("value",1.0);
	allPCcoefs(0,i)=val;
  }

      if(noFileFlag){

  // Reset the counter
  ic=0;
  for (int i=0;i<np;i++){
    
    if(!strcmp(types(i).c_str(),"uncertain")){
	for(int ip=1;ip<=(*outOrder);ip++){
	  char buff[100];
	  sprintf(buff, "cf_%d", ip);
	  string cf_str = buff;
	  Array1D<int> mi(uncdim,0);
	  mi(ic)=ip;
	  allPCcoefs(get_invmindex(mi),i)=params->get_child(i)->attributes()->get_double(cf_str,0.0);
	}

      
      
      ic++;
    }
  
}
      }

      else{
	Array2D<double> uncPCcoefs(npc,uncdim,0.e0);
	read_datafile(uncPCcoefs,inpcfile.c_str());
	for(int ic=0;ic<uncdim;ic++)
	  for(int j=0;j<npc;j++)
	    allPCcoefs(j,uncParamInd(ic))=uncPCcoefs(j,ic);
	}



  return;
}




void readXMLChainInput(RefPtr<XMLElement> xmlTree,MCMC* pmchain, Array1D<double>& chstart, int* nsteps, Array1D<int>& chainParamInd, Array1D<string>& priortype, Array1D<double>& priorparam1, Array1D<double>& priorparam2)
{

  // Select the tree element with the input parameters for this run
  RefPtr<XMLElement> mcmcInput = xmlTree->get_child("MCMC");
  string name = mcmcInput->attributes()->get("name","no_name_specified");
  string method = mcmcInput->attributes()->get("method","no_method_specified");
  string chainparam = mcmcInput->attributes()->get("chainparam","no_chainparam_specified");
  string output = mcmcInput->attributes()->get("output","no_output_specified");
  cout << "Name       : " << name << endl;
  cout << "Method     : " << method << endl;
  cout << "Chainparam : " << chainparam << endl;
  cout << "Output     : " << output << endl;

  RefPtr<XMLElement> method_types = mcmcInput->get_child("method_types");
  RefPtr<XMLElement> pt1=method_types->get_child(method);  
  double gamma=pt1->attributes()->get_double("gamma",1.);
  double eps_cov=pt1->attributes()->get_double("eps_cov",1e-8);
 
  int adaptstart=pt1->attributes()->get_int("adstart",1000);
  int adaptstep=pt1->attributes()->get_int("adstep",10);
  int adaptend=pt1->attributes()->get_int("adstop",1000000);
  (*nsteps)=pt1->attributes()->get_int("nsteps",10000);


  RefPtr<XMLElement> chainparam_types = mcmcInput->get_child("chainparam_types");
  RefPtr<XMLElement> pt2=chainparam_types->get_child(chainparam);  

  Array1D<double> chsig;
  bool noFileFlag;

  if(!strcmp(chainparam.c_str(),"file")) {
    noFileFlag=false;
    /// \todo Put proper checks on the file sizes
    read_datafileVS(chstart, pt2->attributes()->get("chstart","chain_start.dat").c_str());
    read_datafileVS(chsig,   pt2->attributes()->get("chsig","chain_sig.dat").c_str());
  }
  else
    noFileFlag=true;

  RefPtr<XMLElement> params = xmlTree->get_child("ModelParameters");
  int np=params->count_children();
  
  chainParamInd.Clear();
  priortype.Clear();
  priorparam1.Clear();
  priorparam2.Clear();

  for (int i=0;i<np;i++){
    
    string type=params->get_child(i)->attributes()->get("type","no_type_given");
    if(!strcmp(type.c_str(),"infer")){
      chainParamInd.PushBack(i);
      if(noFileFlag){
	double val=params->get_child(i)->attributes()->get_double("value",1.0);
	chstart.PushBack(val);
	chsig.PushBack(params->get_child(i)->attributes()->get_double("sigma",MAX(fabs(0.01*val),0.001)));
      }
      
      priortype.PushBack(params->get_child(i)->attributes()->get("prior","uniform"));
      priorparam1.PushBack(params->get_child(i)->attributes()->get_double("pr1",0.0));
      priorparam2.PushBack(params->get_child(i)->attributes()->get_double("pr2",1000.0));
    }
  }

  int chdim=chainParamInd.XSize();



  RefPtr<XMLElement> output_types = mcmcInput->get_child("output_types");
  RefPtr<XMLElement> pt3=output_types->get_child(output);  
  string filename=pt3->attributes()->get("file","chain.dat");
  int freq_outscreen=pt3->attributes()->get_int("screen",1000);
  int freq_chainfile=pt3->attributes()->get_int("freq",1000);

  pmchain->setChainDim(chdim);
  pmchain->initMethod(method);

  // \todo Need to set checks since the next three lines only useful for method=am
  pmchain->initAdaptSteps(adaptstart,adaptstep,adaptend);
  pmchain->initAMGamma(gamma);
  pmchain->initEpsCov(eps_cov);

  pmchain->setOutputInfo(output,filename,freq_chainfile, freq_outscreen);
  pmchain->initChainPropCovDiag(chsig);
  

  return;
}
