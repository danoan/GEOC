#include "GeneralAdapter.h"

using namespace GEOC::Adapter::GeneralAdapter;

template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
Identity< TIterator, TEstimator, closedCurve>::Identity(MyIterator begin,
                                                        MyIterator end,
                                                        const KSpace& KImage,
                                                        std::vector<EstimationValue>& estimations,
                                                        double h,
                                                        void* data)
{
    AdapterFunctor myFunctor(KImage);
    MyRangeAdapter range(begin,
                         end,
                         myFunctor);

    MyEstimator(range,
                estimations,
                h,
                data);

};

template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
Symmetric< TIterator, TEstimator, closedCurve>::Symmetric(TIterator begin,
                                                          TIterator end,
                                                          const KSpace& KImage,
                                                          std::vector<EstimationValue>& estimations,
                                                          double h,
                                                          void* data)
{
    Curve inverseCurve;
    DIPaCUS::Misc::invertCurve(inverseCurve,
                               KImage,
                               begin,
                               end);

    AdapterFunctor myFunctor(KImage);
    MyRangeAdapter rangeDirectCurve(begin,
                                    end,
                                    myFunctor);

    MyRangeAdapter rangeInverseCurve(inverseCurve.begin(),
                                     inverseCurve.end(),
                                     myFunctor);


    std::vector<EstimationValue> directEstimations;
    std::vector<EstimationValue> inverseEstimations;

    MyEstimator(rangeDirectCurve,
                directEstimations,
                h,
                data);

    MyEstimator(rangeInverseCurve,
                inverseEstimations,
                h,
                data);

    int ip=0;
    int nL = inverseEstimations.size()-1;
    do{
        estimations.push_back( ( directEstimations[ip] - inverseEstimations[nL-ip] )/2.0 );
        ++ip;
    }while(ip<=nL);
};
