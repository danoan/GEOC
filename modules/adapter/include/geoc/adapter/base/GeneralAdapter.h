#ifndef GEOC_ADAPTER_BASE_GENERALADAPTER_H
#define GEOC_ADAPTER_BASE_GENERALADAPTER_H

#include <DGtal/helpers/StdDefs.h>
#include "DIPaCUS/derivates/Misc.h"

#include "geoc/adapter/deducers/EstimatorDeducer.h"
#include "geoc/adapter/deducers/ParameterDeducer.h"

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
            class Identity
            {
            public:
                typedef TIterator MyIterator;

                typedef ParameterDeducer<TEstimator> MyParameterDeducer;
                typedef typename MyParameterDeducer::AdapterFunctor AdapterFunctor;
                typedef typename MyParameterDeducer::EstimationValue EstimationValue;

                typedef EstimatorDeducer<TIterator,TEstimator,AdapterFunctor,EstimationValue,closedCurve> MyEstimator;
                typedef typename MyEstimator::RangeAdapter MyRangeAdapter;
            public:
                Identity(){}

                Identity(MyIterator begin,
                         MyIterator end,
                         const KSpace& KImage,
                         std::vector<EstimationValue>& estimations,
                         double h,
                         void* data);

            };


            template< typename TIterator, template<typename> class TEstimator, bool closedCurve >
            class Symmetric
            {
            public:
                typedef TIterator MyIterator;

                typedef ParameterDeducer<TEstimator> MyParameterDeducer;
                typedef typename MyParameterDeducer::AdapterFunctor AdapterFunctor;
                typedef typename MyParameterDeducer::EstimationValue EstimationValue;

                typedef EstimatorDeducer<TIterator,TEstimator,AdapterFunctor,EstimationValue,closedCurve> MyEstimator;
                typedef typename MyEstimator::RangeAdapter MyRangeAdapter;
            public:
                Symmetric(){}

                Symmetric(MyIterator begin,
                          MyIterator end,
                          const KSpace& KImage,
                          std::vector<EstimationValue>& estimations,
                          double h,
                          void* data);

            };

        }
    }
}

#include "GeneralAdapter.hpp"

#endif //GEOC_ADAPTER_BASE_GENERALADAPTER_H
