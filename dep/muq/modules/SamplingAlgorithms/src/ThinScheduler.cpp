#include "MUQ/SamplingAlgorithms/ThinScheduler.h"

using namespace muq;
using namespace SamplingAlgorithms;

ThinScheduler::ThinScheduler(boost::property_tree::ptree& pt) {
  thinIncr = pt.get("ThinIncrement", 1);
}

bool ThinScheduler::ShouldSave(int step) {
  if (step % thinIncr == 0)
    return true;
  else
    return false;
}
