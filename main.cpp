#include <iostream>
#include <fstream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib> //for exit()

typedef std::vector<double> vd;

double c, dx, dt, alpha, xMax, xMin, tMax, a, A, x_0, pause;
int nxSteps, ntSteps;

void readPars(char*);
void init(vd&, vd&, vd&, vd&);
void firstStep(vd&, vd&, vd&, vd&, vd&, vd&);
void oneStep(vd&, vd&, vd&, vd&, vd&, vd&);
void writeData(vd&, vd&, double);
void writeAnimationScript();

/*==========================================================================*/

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <parameter file>" << std::endl;
        exit(1);
    }

    readPars(argv[1]);

    vd x(nxSteps+1,0);
    vd u0(nxSteps+1,0);
    vd u1(nxSteps+1,0);
    vd r0(nxSteps+1,0);
    vd r1(nxSteps+1,0);
    vd s0(nxSteps+1,0);
    vd s1(nxSteps+1,0);
    double time {0.0};

    init(u0, r0, s0, x);
    firstStep(u0, u1, r0, r1, s0, s1);
    writeData(x, u0, time);

    for(int t=1; t <= ntSteps; t++)
    {
        oneStep(u0, u1, r0, r1, s0, s1);
        time += dt;
        writeData(x, u0, time);
    }
	
	writeAnimationScript();
    return 0;
}

/*==========================================================================*/

void readPars(char* parf)
    {
        std::ifstream parstream(parf);

        if(!parstream)
        {
            std::cerr << "Uh oh, could not open " << parf << " for reading." 
            << std::endl;
            exit(1);
        }

        std::string line;

        while(parstream)
        {
            std::getline(parstream, line);
            if(line.substr(0,3) == "c =") {c = std::stod(line.substr(4));}
            if(line.substr(0,4) == "dt =") {dt = std::stod(line.substr(5));}
            if(line.substr(0,4) == "dx =") {dx = std::stod(line.substr(5));}
            if(line.substr(0,6) == "xMax =") {xMax = std::stod(line.substr(7));}
            if(line.substr(0,6) == "xMin =") {xMin = std::stod(line.substr(7));}
            if(line.substr(0,6) == "tMax =") {tMax = std::stod(line.substr(7));}
            if(line.substr(0,3) == "a =") {a = std::stod(line.substr(4));}
            if(line.substr(0,5) == "x_0 =") {x_0 = std::stod(line.substr(6));}
            if(line.substr(0,3) == "A =") {A = std::stod(line.substr(4));}
            if(line.substr(0,7) == "pause =") {
            	pause = std::stod(line.substr(8));
            }
        }

        alpha = c*dt/dx;
        nxSteps = static_cast<int>((xMax-xMin)/dx);
        ntSteps = static_cast<int>(tMax/dt);
    }

/*==========================================================================*/

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

/*==========================================================================*/

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

/*==========================================================================*/

void oneStep(vd &u0, vd &u1, vd &r0, vd &r1, vd &s0, vd &s1)
{
    for(int i=1; i <= nxSteps; i++)
    {
        s1[i] = s0[i] + alpha*( r0[i] - r0[i-1] );
    }
    s1[0] = s0[0] + alpha*( r0[0] - r0[nxSteps] );

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

/*==========================================================================*/

void writeData(vd &x, vd &u0, double time)
{
	std::ofstream dataf("data.dat", (time == 0) ? std::ios::out : 
						std::ios::app);

	if(!dataf)
	{
		std::cerr << "Uh oh, could not open data.dat for writing" << std::endl;
		exit(1);
	}
	dataf << std::scientific << std::setprecision(14);
	dataf << "#time = " << time << "\n";

	for(int i = 0; i <= nxSteps; i++)
	{
		dataf << x[i] << "\t" << u0[i] << "\n";
	}

	dataf << "\n";
}

/*==========================================================================*/

void writeAnimationScript()
{
	std::ofstream animf("animation.gpi");

	if(!animf)
	{
		std::cerr << "Uh oh, could not open animation.gpi for writing" 
				<< std::endl;
		return;
	}

	animf << "set xrange [" << static_cast<int>(std::floor(xMin)) << ":" 
		<< static_cast<int>(std::ceil(xMax)) << "]\n";
	animf << "set yrange[" << -0.5*A << ":" << A+0.1 << "]\n";
	animf << "do for [i=0:" << ntSteps << "] {\n\ttime = " << dt << "*i\n";
	animf << "\ttitlevar = sprintf(\"time = %f\", time)\n";
	animf << "\tp 'data.dat' every :::i::i w l title titlevar\n";
	animf << "\tpause 0.02\n}\n";    
}