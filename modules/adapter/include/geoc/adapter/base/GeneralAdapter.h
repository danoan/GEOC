#ifndef GEOC_ADAPTER_BASE_GENERALADAPTER_H
#define GEOC_ADAPTER_BASE_GENERALADAPTER_H

#include <DGtal/helpers/StdDefs.h>
#include "geoc/adapter/patch/SCellToPoint.h"

#include "DIPaCUS/derivates/Misc.h"

#include "EstimatorDeducer.h"

namespace GEOC
{
    namespace Adapter
    {
        namespace GeneralAdapter
        {
            typedef DGtal::Z2i::Point Point;
            typedef DGtal::Z2i::Curve Curve;
            typedef DGtal::Z2i::KSpace KSpace;

            template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
            class IdentityRangeCurvature
            {
            public:
                typedef TIterator MyIterator;
                typedef double EstimationValue;

                struct AdapterFunctor
                {
                    typedef Curve::SCell Output;

                    inline
                    Output operator()(const Output& aT) const
                    {
                        return aT;
                    }
                };

                typedef EstimatorDeducer<TIterator,TEstimator,AdapterFunctor,EstimationValue,closedCurve> MyEstimator;
                typedef typename MyEstimator::RangeAdapter MyRangeAdapter;
            public:
                IdentityRangeCurvature(){}

                IdentityRangeCurvature(MyIterator begin,
                                       MyIterator end,
                                       const KSpace& KImage,
                                       std::vector<EstimationValue>& estimations,
                                       double h);

            };


            template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
            class SymmetricCurvature
            {
            public:
                typedef TIterator MyIterator;
                typedef double EstimationValue;

                typedef DGtal::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

                typedef EstimatorDeducer<TIterator,TEstimator,AdapterFunctor,EstimationValue,closedCurve> MyEstimator;
                typedef typename MyEstimator::RangeAdapter MyRangeAdapter;
            public:
                SymmetricCurvature(){}

                SymmetricCurvature(MyIterator begin,
                                   MyIterator end,
                                   const KSpace& KImage,
                                   std::vector<EstimationValue>& estimations,
                                   double h);

            };

            template<typename TIterator, template<typename> class TEstimator, bool closedCurve>
            class SymmetricTangent
            {
            public:
                typedef TIterator MyIterator;
                typedef DGtal::PointVector<2,double> EstimationValue;

                typedef DGtal::Patch::functors::SCellToPoint<KSpace> AdapterFunctor;

                typedef EstimatorDeducer<TIterator,TEstimator,AdapterFunctor,EstimationValue,closedCurve> MyEstimator;
                typedef typename MyEstimator::RangeAdapter MyRangeAdapter;

            public:
                SymmetricTangent(){}

                SymmetricTangent(MyIterator begin,
                                 MyIterator end,
                                 const KSpace& KImage,
                                 std::vector< EstimationValue >& estimationsTangent,
                                 double h);
            };

            template<typename TTangentAdapter>
            class ProjectedLength
            {
            public:
                typedef typename TTangentAdapter::MyIterator MyIterator;
                typedef typename TTangentAdapter::EstimationValue MyTangentValue;

                typedef double EstimationValue;

                typedef TTangentAdapter MyTangentAdapter;

            public:
                ProjectedLength(MyIterator begin,
                                MyIterator end,
                                const KSpace& KImage,
                                std::vector< EstimationValue >& estimations,
                                double h);
            };


            template<typename TTangentAdapter>
            class SinCosLength
            {
            public:
                typedef typename TTangentAdapter::MyIterator MyIterator;
                typedef typename TTangentAdapter::EstimationValue MyTangentValue;

                typedef double EstimationValue;

                typedef TTangentAdapter MyTangentAdapter;

            public:
                SinCosLength(MyIterator begin,
                             MyIterator end,
                             const KSpace &KImage,
                             std::vector<EstimationValue> &estimations,
                             double h);
            };
        }
    }
}

#include "GeneralAdapter.hpp"

#endif //GEOC_ADAPTER_BASE_GENERALADAPTER_H
