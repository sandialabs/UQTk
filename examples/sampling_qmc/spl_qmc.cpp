#include "sampling.hpp"
#include "arrayio.h"

char *getCmdOption(char **begin, char **end, const std::string &option) {

  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    return *itr;

  return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

int main(int argc, char *argv[]) {

  int nspl = 40;
  int ndim = 2;
  char *fout = (char *) "qmc.dat";

  // command-line options
  if( cmdOptionExists(argv, argv+argc,"-d"))         // no. of dimensions
    ndim=atoi(getCmdOption(argv,argv + argc,"-d"));
  if( cmdOptionExists(argv, argv+argc,"-n"))         // no. of samples
    nspl=atoi(getCmdOption(argv,argv + argc,"-n"));
  if( cmdOptionExists(argv, argv+argc,"-f"))         // output filename
    fout=getCmdOption(argv,argv + argc,"-f");

  std::cout<<"No. of dimensions: "<<ndim <<std::endl;
  std::cout<<"No. of samples   : "<<nspl <<std::endl;
  std::cout<<"Output filename  : "<<fout <<std::endl;

  std::vector<double> seq(ndim*nspl);
  Sampling spl(std::string("qmc"),ndim);
  spl.getHaltonSeq ( nspl, ndim, seq );
  write_datafile(seq, nspl, ndim, (char *) "R", fout, (char *)"w");

  return (0);

}
