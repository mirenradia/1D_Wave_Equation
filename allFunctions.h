#ifndef ALL_FUNCTIONS_H
#define ALL_FUNCTIONS_H
#include <vector>

typedef std::vector<double> vd;

//Step takers
void init(vd&, vd&, vd&, vd&);
void firstStep(vd&, vd&, vd&, vd&, vd&, vd&);
void oneStep(vd&, vd&, vd&, vd&, vd&, vd&);

#endif
