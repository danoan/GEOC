#include "DGtal/base/ConstAlias.h"


namespace DGtal
{
    namespace Patch
    {
        namespace functors
        {
            template <typename KSpace>
            class SCellToPoint
            {
            public:
                typedef typename KSpace::Point Output;
                typedef typename KSpace::SCell Input;

            private:
                /**
                  * Aliasing pointer on the Khalimsky space.
                 */
                const KSpace* myK;

            public:

                /**
                  * Default constructor.
                 */
                SCellToPoint() : myK(NULL) { }
                /**
                  * Constructor.
                  * @param aK a Khalimsky space
                 */
                SCellToPoint( ConstAlias<KSpace> aK ) : myK(&aK) { }

                /**
                 * Copy constructor.
                 * @param other any SCellToPoint functor
                 */
                SCellToPoint(const SCellToPoint& other)
                        : myK(other.myK) { }

                /**
                 * Assignment.
                 *
                 * @param other the object to copy.
                 * @return a reference on 'this'.
                 */
                SCellToPoint& operator= ( const SCellToPoint & other )
                {
                    if (this != &other)
                    {
                        myK = other.myK;
                    }
                    return *this;
                }

                /**
                 * Returns a point (with integer coordinates)
                 * from a scell (with khalimsky coordinates)
                 * @param aSCell a scell
                 * @return the corresponding point.
                 */
                Output operator()(const Input& aSCell) const
                {
                    ASSERT( myK );
                    Input s = aSCell;
                    while ( myK->sDim(s) > 0 )
                    {
                        Input tmp( myK->sIndirectIncident( s, *myK->sDirs( s ) ) );
                        ASSERT( myK->sDim(tmp) < myK->sDim(s) );
                        s = tmp;
                    }
                    return Output( myK->sCoords(s) );
                }

            }; // end of class SCellToPoint
        }
    }
}