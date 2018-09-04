#include "GridCurve.h"

using namespace GEOC::Adapter::GridCurve;

template<typename TEstimator>
SymmetricCurvature<TEstimator>::SymmetricCurvature(Curve::ConstIterator begin,
                                                   Curve::ConstIterator end,
                                                   const KSpace& KImage,
                                                   std::vector<double>& estimations,
                                                   bool closedCurve)
{
    Curve negativeCurve;
    DIPaCUS::Misc::InvertCurve(KImage,
                               begin,
                               end,
                               negativeCurve);

    AdapterFunctor myFunctor(KImage);
    SCellToPointRangeAdapter rangePositiveCurve(begin,
                                                end,
                                                myFunctor);

    SCellToPointRangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                                negativeCurve.end(),
                                                myFunctor);


    std::vector<double> positiveEstimations;
    std::vector<double> negativeEstimations;

    if(closedCurve) {
        TEstimator(rangePositiveCurve.c(),
                  rangePositiveCurve.c(),
                  positiveEstimations);

        TEstimator(rangeNegativeCurve.c(),
                  rangeNegativeCurve.c(),
                  negativeEstimations);
    }else{
        TEstimator(rangePositiveCurve.begin(),
                  rangePositiveCurve.end(),
                  positiveEstimations);

        TEstimator(rangeNegativeCurve.begin(),
                  rangeNegativeCurve.end(),
                  negativeEstimations);
    }

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( ( fabs(positiveEstimations[ip]) + fabs(negativeEstimations[nL-ip]) )/2.0 );
        ++ip;
    }while(ip<=nL);
};

template<typename TEstimator>
SymmetricTangent<TEstimator>::SymmetricTangent(Curve::ConstIterator begin,
                                               Curve::ConstIterator end,
                                               const KSpace& KImage,
                                               std::vector< TangentVector >& estimationsTangent,
                                               bool closedCurve)
{
    AdapterFunctor myFunctor(KImage);
    SCellToPointRangeAdapter rangePositiveCurve(begin,
                                                end,
                                                myFunctor);

    Curve negativeCurve;
    DIPaCUS::Misc::InvertCurve(KImage,begin,end,negativeCurve);

    SCellToPointRangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                                negativeCurve.end(),
                                                myFunctor);

    std::vector<TangentVector> positiveEstimations;
    std::vector<TangentVector> negativeEstimations;

    if(closedCurve) {

        TEstimator(rangePositiveCurve.c(),
                   rangePositiveCurve.c(),
                   positiveEstimations);

        TEstimator(rangeNegativeCurve.c(),
                   rangeNegativeCurve.c(),
                   negativeEstimations);
    }else{
        TEstimator(rangePositiveCurve.begin(),
                   rangePositiveCurve.end(),
                   positiveEstimations);

        TEstimator(rangeNegativeCurve.begin(),
                   rangeNegativeCurve.end(),
                   negativeEstimations);
    }


    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimationsTangent.push_back( (positiveEstimations[ip]-negativeEstimations[nL-ip]).getNormalized() );
        ++ip;
    }while(ip<=nL);
}
