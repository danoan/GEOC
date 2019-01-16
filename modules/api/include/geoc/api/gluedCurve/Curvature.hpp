#ifndef GEOC_API_GLUEDCURVE_CURVATURE_H
#define GEOC_API_GLUEDCURVE_CURVATURE_H

#include <DGtal/helpers/StdDefs.h>

#include <geoc/estimator/adaptable/Curvature.h>
#include <geoc/adapter/GluedCurve.h>

namespace GEOC
{
    namespace API
    {
        namespace GluedCurve
        {
            class Curvature
            {
            public:

                template<typename T>
                using ALG_MDCA = GEOC::Estimator::Standard::MDCACurvature<T>;

                template<typename T>
                using ALG_II = GEOC::Estimator::Standard::IICurvature<T>;

                template<typename T>
                using ALG_LMDCA = GEOC::Estimator::Alternative::LMDCACurvature<T>;

                template<typename T>
                using ALG_HMDCA = GEOC::Estimator::Alternative::HMDCACurvature<T>;

            public:
                typedef DGtal::Z2i::KSpace KSpace;
                typedef GEOC::Adapter::GluedCurve::IteratorType CurveIterator;

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
                    GEOC::Adapter::GluedCurve::Symmetric< TAlgorithm,false>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                static void symmetricClosed(const KSpace& KImage,
                                            CurveIterator begin,
                                            CurveIterator end,
                                            EstimationsVector& ev,
                                            double h=1.0)
                {
                    GEOC::Adapter::GluedCurve::Symmetric< TAlgorithm,true>(begin,end,KImage,ev,h);
                }

            };
        }
    }
}

#endif //GEOC_API_GLUEDCURVE_CURVATURE_H
