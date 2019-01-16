#include "GeneralAdapter.h"

using namespace GEOC::Adapter::GeneralAdapter;

template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
Identity< TIterator, TEstimator, closedCurve>::Identity(MyIterator begin,
                                                        MyIterator end,
                                                        const KSpace& KImage,
                                                        std::vector<EstimationValue>& estimations,
                                                        double h)
{
    AdapterFunctor myFunctor;
    MyRangeAdapter range(begin,
                         end,
                         myFunctor);

    MyEstimator(range,
                estimations,
                h);

};

template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
Symmetric< TIterator, TEstimator, closedCurve>::Symmetric(TIterator begin,
                                                          TIterator end,
                                                          const KSpace& KImage,
                                                          std::vector<EstimationValue>& estimations,
                                                          double h)
{
    Curve negativeCurve;
    DIPaCUS::Misc::InvertCurve<Curve::ConstIterator>(KImage,
                                                     begin,
                                                     end,
                                                     negativeCurve);

    AdapterFunctor myFunctor(KImage);
    MyRangeAdapter rangePositiveCurve(begin,
                                             end,
                                             myFunctor);

    MyRangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                             negativeCurve.end(),
                                             myFunctor);


    std::vector<EstimationValue> positiveEstimations;
    std::vector<EstimationValue> negativeEstimations;

    MyEstimator(rangePositiveCurve,
                positiveEstimations,
                h);

    MyEstimator(rangeNegativeCurve,
                negativeEstimations,
                h);

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( ( fabs(positiveEstimations[ip]) + fabs(negativeEstimations[nL-ip]) )/2.0 );
        ++ip;
    }while(ip<=nL);
};
