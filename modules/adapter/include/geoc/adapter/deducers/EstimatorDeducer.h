#ifndef GEOC_ADAPTER_ESTIMATORDEDUCER_H
#define GEOC_ADAPTER_ESTIMATORDEDUCER_H

#include <DGtal/helpers/StdDefs.h>

namespace GEOC
{
    namespace Adapter
    {
        template< typename TIterator, template<typename> class TEstimator, typename TAdapterFunctor, typename TValue, bool closedCurve >
        class EstimatorDeducer
        {
        public:
            typedef DGtal::Z2i::KSpace KSpace;
            typedef TAdapterFunctor AdapterFunctor;

            typedef DGtal::ConstRangeAdapter< TIterator,
                    AdapterFunctor,
                    typename AdapterFunctor::Output > RangeAdapter;

        public:
            EstimatorDeducer(RangeAdapter range,
                             std::vector<TValue>& estimations,
                             double h);
        };

        template< typename TIterator, template<typename> class TEstimator, typename TAdapterFunctor, typename TValue>
        class EstimatorDeducer<TIterator,TEstimator,TAdapterFunctor,TValue,false>
        {
        public:
            typedef DGtal::Z2i::KSpace KSpace;
            typedef TAdapterFunctor AdapterFunctor;

            typedef DGtal::ConstRangeAdapter< TIterator,
                    AdapterFunctor,
                    typename AdapterFunctor::Output > RangeAdapter;

            typedef typename RangeAdapter::ConstIterator RangeIterator;
            typedef TEstimator<RangeIterator> MyEstimator;

        public:
            EstimatorDeducer(RangeAdapter range,
                             std::vector<TValue>& estimations,
                             double h)
            {
                MyEstimator(range.begin(),range.end(),estimations,h);
            }
        };

        template< typename TIterator, template<typename> class TEstimator, typename TAdapterFunctor, typename TValue>
        class EstimatorDeducer<TIterator,TEstimator,TAdapterFunctor,TValue,true>
        {
        public:
            typedef DGtal::Z2i::KSpace KSpace;
            typedef TAdapterFunctor AdapterFunctor;

            typedef DGtal::ConstRangeAdapter< TIterator,
                    AdapterFunctor,
                    typename AdapterFunctor::Output > RangeAdapter;

            typedef typename RangeAdapter::ConstCirculator RangeCirculator;
            typedef TEstimator<RangeCirculator> MyEstimator;

            EstimatorDeducer(RangeAdapter range,
                             std::vector<TValue>& estimations,
                             double h)
            {
                MyEstimator(range.c(),range.c(),estimations,h);
            }
        public:
        };
    }
}

#endif //GEOC_ADAPTER_ESTIMATORDEDUCER_H
