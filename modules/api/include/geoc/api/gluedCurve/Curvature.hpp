#ifndef GEOC_GLUEDCURVECURVATURE_H
#define GEOC_GLUEDCURVECURVATURE_H

#include <DGtal/helpers/StdDefs.h>

#include <geoc/estimator/standard/Curvature.h>
#include <geoc/estimator/alternative/Curvature.h>
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
                using ALG_LDCA = GEOC::Estimator::Alternative::DCALambdaCurvature<T>;

                template<typename T>
                using ALG_HDCA = GEOC::Estimator::Alternative::HDCACurvature<T>;

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
                    GEOC::Adapter::GluedCurve::SymmetricCurvature< TAlgorithm,false>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                static void symmetricClosed(const KSpace& KImage,
                                            CurveIterator begin,
                                            CurveIterator end,
                                            EstimationsVector& ev,
                                            double h=1.0)
                {
                    GEOC::Adapter::GluedCurve::SymmetricCurvature< TAlgorithm,true>(begin,end,KImage,ev,h);
                }

            };
        }
    }
}

#endif //GEOC_GRIDCURVECURVATURE_H
