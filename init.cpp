#include "constants.h"
#include <vector>
#include <cmath>

typedef std::vector<double> vd;

void init(vd &u0, vd &r0, vd &s0, vd &x)
{
    for(int i=0; i <= nxSteps; i++)
    {
        x[i] = xMin + i*dx;
        u0[i] = A*std::exp(-(x[i] - x_0)*(x[i] - x_0)/(2*a*a));
    }

    for(int i=0; i < nxSteps; i++)
    {
        r0[i] = (c/dx)*(u0[i+1] - u0[i]);
    }
    r0[nxSteps] = (c/dx)*(u0[0] - u0[nxSteps]);

    s0 = vd(nxSteps+1,0);
}
