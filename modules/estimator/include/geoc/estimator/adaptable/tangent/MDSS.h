#ifndef GEOC_ESTIMATORS_CURVATURE_MDSS_H
#define GEOC_ESTIMATORS_CURVATURE_MDSS_H

#include <DGtal/geometry/curves/ArithmeticalDSSComputer.h>
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include <DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h>

namespace GEOC
{
    namespace Estimator
    {
        namespace Standard
        {
            template<typename IteratorType>
            struct MDSSTangent
            {
                typedef DGtal::ArithmeticalDSSComputer< IteratorType, int, 4 > SegmentComputer;
                typedef DGtal::TangentFromDSSEstimator<SegmentComputer> SCEstimator;

                typedef DGtal::PointVector<2,double> TangentVector;

                MDSSTangent(IteratorType itb,
                            IteratorType ite,
                            std::vector< TangentVector >& estimations,
                            double h,
                            void* data)
                {
                    SegmentComputer sc;
                    SCEstimator f;

                    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDSSTangentEstimator(sc,f);

                    MCMDSSTangentEstimator.init(h,itb,ite);
                    MCMDSSTangentEstimator.eval(itb,ite,std::back_inserter(estimations));
                }
            };
        }
    }
}

#endif //GEOC_ESTIMATORS_CURVATURE_MDSS_H
