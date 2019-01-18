#ifndef GEOC_ADAPTERFUNCTORSOLVER_H
#define GEOC_ADAPTERFUNCTORSOLVER_H

#include <DGtal/helpers/StdDefs.h>
#include <geoc/adapter/patch/SCellToPoint.h>

#include <geoc/estimator/adaptable/Curvature.h>
#include <geoc/estimator/adaptable/Tangent.h>


namespace GEOC
{
    namespace Adapter
    {
        template<typename T>
        using MDCA = GEOC::Estimator::Standard::MDCACurvature<T>;

        template<typename T>
        using II = GEOC::Estimator::Standard::IICurvature<T>;
        
        template<typename T>
        using HMDCA = GEOC::Estimator::Alternative::HMDCACurvature<T>;

        template<typename T>
        using LMDCA = GEOC::Estimator::Alternative::LMDCACurvature<T>;


        template<typename T>
        using MDSS = GEOC::Estimator::Standard::MDSSTangent<T>;

        template<typename T>
        using LMDSS = GEOC::Estimator::Alternative::LMDSSTangent<T>;
        
        
        template< template<typename> class TEstimator>
        struct ParameterDeducer
        {
            typedef double EstimationValue;
            struct AdapterFunctor
            {
                typedef DGtal::Z2i::KSpace KSpace;
                typedef DGtal::Z2i::Curve::SCell Output;

                AdapterFunctor(const KSpace& KImage){}

                inline
                Output operator()(const Output& aT) const
                {
                    return aT;
                }
            };
        };

        template<>
        struct ParameterDeducer<MDCA>
        {
            typedef double EstimationValue;
            typedef DGtal::functors::SCellToIncidentPoints<DGtal::Z2i::KSpace> AdapterFunctor;
        };

        template<>
        struct ParameterDeducer<HMDCA>
        {
            typedef double EstimationValue;
            typedef DGtal::functors::SCellToIncidentPoints<DGtal::Z2i::KSpace> AdapterFunctor;
        };

        template<>
        struct ParameterDeducer<LMDCA>
        {
            typedef double EstimationValue;
            typedef DGtal::functors::SCellToIncidentPoints<DGtal::Z2i::KSpace> AdapterFunctor;
        };        

        template<>
        struct ParameterDeducer<MDSS>
        {
            typedef DGtal::PointVector<2,double> EstimationValue;
            typedef DGtal::Patch::functors::SCellToPoint<DGtal::Z2i::KSpace> AdapterFunctor;
        };

        template<>
        struct ParameterDeducer<LMDSS>
        {
            typedef DGtal::PointVector<2,double> EstimationValue;
            typedef DGtal::Patch::functors::SCellToPoint<DGtal::Z2i::KSpace> AdapterFunctor;
        };        
    }
}

#endif //GEOC_ADAPTERFUNCTORSOLVER_H
