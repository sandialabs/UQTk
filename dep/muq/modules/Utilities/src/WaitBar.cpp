#include "MUQ/Utilities/WaitBar.h"
#include <stdio.h>
#include <iostream>

using namespace muq::Utilities;

WaitBar::WaitBar(double minValIn, double maxValIn) : minVal(minValIn), maxVal(maxValIn)
{
  Update(minVal);
}

void WaitBar::Complete()
{
  std::cout << std::endl;
}

void WaitBar::Update(double newVal)
{
  double perc = (newVal-minVal)/(maxVal-minVal);
  int pos = double(barWidth) * perc ;

  printf("[");
  for(int i=0; i<pos; ++i)
    printf("=");

  printf(">");
  
  for(int i=pos+1; i<barWidth; ++i)
    printf(" ");

  printf("] %3.1f%% \r", (perc*100.0));

  std::cout << std::flush;
}
