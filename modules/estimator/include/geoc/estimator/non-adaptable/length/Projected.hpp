#ifndef GEOC_ESTIMATOR_LENGTH_PROJECTED_H
#define GEOC_ESTIMATOR_LENGTH_PROJECTED_H

#include <DGtal/helpers/StdDefs.h>

#include <vector>
#include <cmath>

namespace GEOC
{
    namespace Estimator
    {
        namespace Standard
        {
            template<typename TTangentAdapter>
            class ProjectedLength
            {
            public:
                typedef DGtal::Z2i::Point Point;
                typedef DGtal::Z2i::KSpace KSpace;

                typedef typename TTangentAdapter::MyIterator MyIterator;
                typedef typename TTangentAdapter::EstimationValue MyTangentValue;

                typedef double EstimationValue;
                typedef TTangentAdapter MyTangentAdapter;

            public:

                ProjectedLength(MyIterator begin,
                                MyIterator end,
                                const KSpace& KImage,
                                std::vector< EstimationValue >& estimations,
                                double h)
                {
                    std::vector< MyTangentValue > tangentEstimations;
                    MyTangentAdapter(begin,end,KImage,tangentEstimations,h);

                    Point pTarget,pSource,scellVector;
                    auto it = begin;
                    int i = 0;
                    do
                    {
                        pTarget = KImage.sCoords( KImage.sDirectIncident(*it,*KImage.sDirs(*it)) );
                        pSource = KImage.sCoords( KImage.sIndirectIncident(*it,*KImage.sDirs(*it)) );

                        scellVector = pTarget-pSource;

                        //Scale by h because the tangent is normalized
                        estimations.push_back( h*fabs( ( tangentEstimations[i].dot(scellVector) ) ) );

                        ++i;
                        ++it;
                    }while(it!=end);
                }
            };
        }
    }
}

#endif //GEOC_ESTIMATOR_LENGTH_PROJECTED_H
