#include <iostream>
#include <ios>
#include <iomanip>
#include "constants.h"
#include "allFunctions.h"

int main()
{
    typedef std::vector<double> vd;

    vd x(nxSteps+1,0);
    vd u0(nxSteps+1,0);
    vd u1(nxSteps+1,0);
    vd r0(nxSteps+1,0);
    vd r1(nxSteps+1,0);
    vd s0(nxSteps+1,0);
    vd s1(nxSteps+1,0);

    init(u0,r0,s0,x);
    firstStep(u0,u1,r0,r1,s0,s1);

    for(int t=1; t < ntSteps; t++)
    {
        oneStep(u0, u1, r0, r1, s0, s1);
    }

    for(int i=0; i <= nxSteps; i++)
    {
        std::cout << std::setprecision(3) << std::fixed;
        std::cout << x[i] << "\t";
        std::cout << std::setprecision(16) << std::scientific;
        std::cout << u0[i] << "\t" << s0[i] << "\t" <<r0[i] << "\n";
    }

    return 0;
}
