using namespace GEOC::Adapter::GluedCurve;

template<typename TEstimator>
SymmetricCurvature<TEstimator>::SymmetricCurvature(GCurve::GluedLinelsIterator begin,
                                                   GCurve::GluedLinelsIterator end,
                                                   const KSpace& KImage,
                                                   std::vector<double>& estimations)
{
    AdapterFunctor myFunctor(KImage);
    GluedRangeAdapter rangePositiveCurve(begin,
                                         end,
                                         myFunctor);

    Curve negativeCurve;
    DIPaCUS::Misc::InvertCurve<GCurve::GluedLinelsIterator>(KImage,
                                                            begin,
                                                            end,
                                                            negativeCurve);


    CurveRangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                         negativeCurve.end(),
                                         myFunctor);

    std::vector<double> positiveEstimations;
    std::vector<double> negativeEstimations;

    TEstimator(rangePositiveCurve.begin(),
               rangePositiveCurve.end(),
               positiveEstimations);

    TEstimator(rangeNegativeCurve.begin(),
               rangeNegativeCurve.end(),
               negativeEstimations);


    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( ( fabs(positiveEstimations[ip]) + fabs(negativeEstimations[nL-ip]) )/2.0 );
        ++ip;
    }while(ip<=nL);
}

template<typename TEstimator>
SymmetricTangent<TEstimator>::SymmetricTangent(GCurve::GluedLinelsIterator begin,
                                               GCurve::GluedLinelsIterator end,
                                               const KSpace& KImage,
                                               std::vector< TangentVector >& estimations)
{
    AdapterFunctor myFunctor(KImage);
    GluedRangeAdapter rangePositiveCurve(begin,
                                         end,
                                         myFunctor);

    Curve negativeCurve;
    DIPaCUS::Misc::InvertCurve<GCurve::GluedLinelsIterator>(KImage,
                                                            begin,
                                                            end,
                                                            negativeCurve);


    CurveRangeAdapter rangeNegativeCurve(negativeCurve.begin(),
                                         negativeCurve.end(),
                                         myFunctor);

    std::vector<TangentVector> positiveEstimations;
    TEstimator(rangePositiveCurve.begin(),
               rangePositiveCurve.end(),
               positiveEstimations);

    std::vector<TangentVector> negativeEstimations;
    TEstimator(rangeNegativeCurve.begin(),
               rangeNegativeCurve.end(),
               negativeEstimations);

    int ip=0;
    int nL = negativeEstimations.size()-1;
    do{
        estimations.push_back( (positiveEstimations[ip]-negativeEstimations[nL-ip]).getNormalized() );
        ++ip;
    }while(ip<=nL);

}


template<typename TTangentAdapter>
ProjectedLength<TTangentAdapter>::ProjectedLength(GCurve::GluedLinelsIterator begin,
                                                  GCurve::GluedLinelsIterator end,
                                                  const KSpace &KImage,
                                                  std::vector<double> &estimations)
{
    TangentVector tangentEstimations;
    TTangentAdapter(begin,end,KImage,tangentEstimations);

    Point pTarget,pSource,scellVector;
    for(auto it=tangentEstimations.begin();it!=tangentEstimations.end();++it)
    {
        pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
        pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );

        scellVector = pTarget-pSource;
        estimations.push_back( fabs( tangentEstimations[i].dot(scellVector) ) );
    }

};


template<typename TTangentAdapter>
SinCosLength<TTangentAdapter>::SinCosLength(GCurve::GluedLinelsIterator begin,
                                            GCurve::GluedLinelsIterator end,
                                            const KSpace &KImage,
                                            std::vector<double> &estimations)
{
    TangentVector tangentEstimations;
    TTangentAdapter(begin,end,KImage,tangentEstimations);

    for(auto it=tangentEstimations.begin();it!=tangentEstimations.end();++it)
    {
        //1.0/(cos+sin) Length estimation
        estimations.push_back( 1.0/( fabs( (*it)[0] ) + fabs( (*it)[1] ) ) );
    }

};