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
    Curve inverseCurve;
    DIPaCUS::Misc::InvertCurve<Curve::ConstIterator>(KImage,
                                                     begin,
                                                     end,
                                                     inverseCurve);

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
                h);

    MyEstimator(rangeInverseCurve,
                inverseEstimations,
                h);

    int ip=0;
    int nL = inverseEstimations.size()-1;
    do{
        estimations.push_back( ( directEstimations[ip] + inverseEstimations[nL-ip] )/2.0 );
        ++ip;
    }while(ip<=nL);
};
