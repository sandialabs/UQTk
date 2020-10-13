#ifndef WAITBAR_H
#define WAITBAR_H


namespace muq{
  namespace Utilities{


    class WaitBar{

    public:
      WaitBar(double minValIn, double maxValIn);
      
      void Complete();

      void Update(double newVal);

    private:
      double minVal, maxVal;

      const int barWidth = 60;
    };

  } // namespace Utilities
}// namespace muq



#endif // WAITBAR_H
