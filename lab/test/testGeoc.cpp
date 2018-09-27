#include <TestConvergence.h>

using namespace GEOC::Test;

namespace GEOC
{
    namespace Test
    {
        bool verbose = true;
        bool visualOutput = false;
    }
}

int main()
{
    TestConvergence(5,1.0);
    std::cout << std::endl << std::endl;

    TestConvergence(5,0.5);
    std::cout << std::endl << std::endl;

    TestConvergence(5,0.25);
    std::cout << std::endl << std::endl;

    TestConvergence(5,0.1);
    std::cout << std::endl << std::endl;

    TestConvergence(5,0.05);
    return 0;
}

