#ifndef GEOC_GRIDCURVECURVATURE_H
#define GEOC_GRIDCURVECURVATURE_H

#include <DGtal/helpers/StdDefs.h>

#include <geoc/estimator/standard/Curvature.h>
#include <geoc/adapter/GridCurve.h>
#include <geoc/estimator/alternative/Curvature.h>

namespace GEOC
{
    namespace API
    {
        namespace GridCurve
        {
            class Curvature
            {
            public:

                template<typename T>
                using ALG_MDCA = GEOC::Estimator::Standard::MDCACurvature<T>;

                template<typename T>
                using ALG_II = GEOC::Estimator::Standard::IICurvature<T>;

                template<typename T>
                using ALG_LDCA = GEOC::Estimator::Alternative::DCALambdaCurvature<T>;

                template<typename T>
                using ALG_HDCA = GEOC::Estimator::Alternative::HDCACurvature<T>;

            public:
                typedef DGtal::Z2i::Curve Curve;
                typedef Curve::ConstIterator CurveIterator;

                typedef std::vector<double> EstimationsVector;

            private:
                Curvature();
                Curvature(const Curvature&)=delete;
                Curvature& operator=(const Curvature&)=delete;

            public:

                template< template<typename> class TAlgorithm>
                static void symmetricOpen(const KSpace& KImage,
                                          CurveIterator begin,
                                          CurveIterator end,
                                          EstimationsVector& ev,
                                          double h=1.0)
                {
                    GEOC::Adapter::GridCurve::SymmetricCurvature< TAlgorithm,false>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                static void symmetricClosed(const KSpace& KImage,
                                            CurveIterator begin,
                                            CurveIterator end,
                                            EstimationsVector& ev,
                                            double h=1.0)
                {
                    GEOC::Adapter::GridCurve::SymmetricCurvature< TAlgorithm,true>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                static void identityOpen(const KSpace& KImage,
                                           CurveIterator begin,
                                           CurveIterator end,
                                           EstimationsVector& ev,
                                           double h=1.0)
                {
                    GEOC::Adapter::GridCurve::IdentityRangeCurvature< TAlgorithm,false>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                static void identityClosed(const KSpace& KImage,
                                            CurveIterator begin,
                                            CurveIterator end,
                                            EstimationsVector& ev,
                                            double h=1.0)
                {
                    GEOC::Adapter::GridCurve::IdentityRangeCurvature< TAlgorithm,true>(begin,end,KImage,ev,h);
                }

            };
        }
    }
}

#endif //GEOC_GRIDCURVECURVATURE_H
