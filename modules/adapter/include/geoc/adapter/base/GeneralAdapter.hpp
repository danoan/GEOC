#include "GeneralAdapter.h"

using namespace GEOC::Adapter::GeneralAdapter;

template< typename TIterator, template<typename> typename TEstimator, bool closedCurve >
SymmetricCurvature< TIterator, TEstimator, closedCurve>::SymmetricCurvature(TIterator begin,
                                                                            TIterator end,
                                                                            const KSpace& KImage,
                                                                            std::vector<EstimationValue>& estimations)
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
                rangeNegativeCurve,
                positiveEstimations,
                negativeEstimations);

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( ( fabs(positiveEstimations[ip]) + fabs(negativeEstimations[nL-ip]) )/2.0 );
        ++ip;
    }while(ip<=nL);
};

template<typename TIterator, template<typename> typename TEstimator, bool closedCurve>
SymmetricTangent<TIterator, TEstimator, closedCurve>::SymmetricTangent(TIterator begin,
                                                                       TIterator end,
                                                                       const KSpace& KImage,
                                                                       std::vector< EstimationValue >& estimationsTangent)
{
    AdapterFunctor myFunctor(KImage);
    MyRangeAdapter rangePositiveCurve(begin,
                                      end,
                                      myFunctor);

    Curve negativeCurve;
    DIPaCUS::Misc::InvertCurve<Curve::ConstIterator>(KImage,
                                                     begin,
                                                     end,
                                                     negativeCurve);

    MyRangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                      negativeCurve.end(),
                                      myFunctor);

    std::vector<EstimationValue> positiveEstimations;
    std::vector<EstimationValue> negativeEstimations;


    MyEstimator(rangePositiveCurve,
                rangeNegativeCurve,
                positiveEstimations,
                negativeEstimations);

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimationsTangent.push_back( (positiveEstimations[ip]-negativeEstimations[nL-ip]).getNormalized() );
        ++ip;
    }while(ip<=nL);
}

template<typename TTangentAdapter>
ProjectedLength<TTangentAdapter>::ProjectedLength(MyIterator begin,
                                                  MyIterator end,
                                                  const KSpace &KImage,
                                                  std::vector<EstimationValue> &estimations)
{
    std::vector< MyTangentValue > tangentEstimations;
    MyTangentAdapter(begin,end,KImage,tangentEstimations);

    Point pTarget,pSource,scellVector;
    auto it = begin;
    int i = 0;
    do
    {
        pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
        pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );

        scellVector = pTarget-pSource;
        estimations.push_back( fabs( ( tangentEstimations[i].dot(scellVector) ) ) );

        ++i;
        ++it;
    }while(it!=end);

};


template<typename TTangentAdapter>
SinCosLength<TTangentAdapter>::SinCosLength(MyIterator begin,
                                            MyIterator end,
                                            const KSpace &KImage,
                                            std::vector<EstimationValue> &estimations)
{
    std::vector< MyTangentValue > tangentEstimations;
    MyTangentAdapter(begin,end,KImage,tangentEstimations);

    auto it = begin;
    int i = 0;
    do
    {
        //1.0/(cos+sin) Length estimation
        estimations.push_back( 1.0/( fabs( tangentEstimations[i][0] ) + fabs( tangentEstimations[i][1] ) ) );
        ++i;
        ++it;
    }while(it!=end);

};
