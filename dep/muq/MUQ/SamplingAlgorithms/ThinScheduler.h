#ifndef THINSCHEDULER_H_
#define THINSCHEDULER_H_

#include "MUQ/SamplingAlgorithms/SaveSchedulerBase.h"

#include <boost/property_tree/ptree.hpp>

namespace muq {
  namespace SamplingAlgorithms {

    class ThinScheduler : public SaveSchedulerBase {

    public:
      ThinScheduler(boost::property_tree::ptree& pt);

      virtual ~ThinScheduler() = default;

      virtual bool ShouldSave(int step) override;

    private:
      int thinIncr;

    };

  } // namespace SamplingAlgoirthms
} // namespace muq

#endif
