#include "constants.h"
#include <vector>

typedef std::vector<double> vd;

void firstStep(vd &u0, vd &u1, vd &r0, vd &r1, vd &s0, vd &s1)
{
    //the 0.5 is because this first step is non-centred
    for(int i=1; i <= nxSteps; i++)
    {
        s1[i] = s0[i] + 0.5*alpha*( r0[i] - r0[i-1] );
    }
    s1[0] = s0[0] + 0.5*alpha*( r0[0] - r0[nxSteps] );

    for(int i=0; i < nxSteps; i++)
    {
        r1[i] = r0[i] + alpha*( s1[i+1] - s1[i] );
    }
    r1[nxSteps] = r0[nxSteps] + alpha*( s1[0] - s1[nxSteps] );

    for(int i=0; i <= nxSteps; i++)
    {
        u1[i] = u0[i] + dt*s1[i];
    }

    u0.swap(u1);
    r0.swap(r1);
    s0.swap(s1);
}
