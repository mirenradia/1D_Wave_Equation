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
int nxSteps, ntSteps, n_icn;
std::string alg;

void readPars(char*);
void init(vd&, vd&, vd&, vd&);
void init_icn(vd&, vd&, vd&, vd&);
void firstStep_leapfrog(vd&, vd&, vd&, vd&, vd&, vd&);
void oneStep_leapfrog(vd&, vd&, vd&, vd&, vd&, vd&);
void oneIter_icn(vd&, vd&, vd&, vd&, vd&, vd&);
void averageIter_icn(vd&, vd&, vd&, vd&);
void calc_u_icn(vd&, vd&, vd&, vd&, vd&, vd&);
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

    //instantiate our vectors to hold solution and first derivatives
    //the 0 quantities are the values at the current time step
    //the 1 quantities are the values at the next time step
    vd x(nxSteps+1,0);
    vd u0(nxSteps+1,0);
    vd u1(nxSteps+1,0);
    vd r0(nxSteps+1,0);
    vd r1(nxSteps+1,0);
    vd s0(nxSteps+1,0);
    vd s1(nxSteps+1,0);
    double time {0.0};

    
    if(alg == "leapfrog")
    { 
        init(u0, r0, s0, x);
        writeData(x, u0, time);
        firstStep_leapfrog(u0, u1, r0, r1, s0, s1);
        time += dt;
        writeData(x, u0, time);

        for(int t=1; t <= ntSteps; t++)
        {
            oneStep_leapfrog(u0, u1, r0, r1, s0, s1);
            time += dt;
            writeData(x, u0, time);
        }
	}
    else if(alg == "icn")
    {
        init_icn(u0, r0, s0, x);
        writeData(x, u0, time);
        
        //instantiate 2 more vectors to hold intermediate values in algorithm
        //these will hold the "1/2" values in the ICN iterations
        vd r01(nxSteps+1,0);
        vd s01(nxSteps+1,0);

        for(int t=0; t <= ntSteps; t++)
        {
            for(int i_icn=0; i_icn <= n_icn; i_icn++)
            {
                if(i_icn == 0)
                {
                    oneIter_icn(r0, r0, r01, s0, s0, s01);
                }
                else if(i_icn != n_icn)
                {
                    averageIter_icn(r0, r01, s0, s01);
                    oneIter_icn(r0, r01, r1, s0, s01, s1);
                    
                    //rotate vectors as this is not the final iteration
                    r1.swap(r01);
                    s1.swap(s01);
                }
                else
                {
                    averageIter_icn(r0, r01, s0, s01);
                    oneIter_icn(r0, r01, r1, s0, s01, s1);
                }
            }
            calc_u_icn(u0, u1, s0, s1, r0, r1);
            time += dt;
            writeData(x, u0, time);
        }
    }
    else
    {
        std::cerr << "In " << argv[1] << " there must be a line with either:\n";
        std::cerr << "alg = icn\nalg = leapfrog" << std::endl;
        return 1;
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
            if(line.substr(0,5) == "alg =") {alg = line.substr(6);}
            if(line.substr(0,7) == "n_icn =") {n_icn = std::stoi(line.substr(8));}
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
    
    //boundary term
    r0[nxSteps] = (c/dx)*(u0[0] - u0[nxSteps]);

    s0 = vd(nxSteps+1,0);
}

/*==========================================================================*/

void init_icn(vd &u0, vd &r0, vd &s0, vd &x)
{
    for(int i=0; i <= nxSteps; i++)
    {
        x[i] = xMin + i*dx;
        u0[i] = A*std::exp(-(x[i] - x_0)*(x[i] - x_0)/(2*a*a));
    }

    for(int i=1; i < nxSteps; i++)
    {
        r0[i] = 0.5*(c/dx)*(u0[i+1] - u0[i-1]);
    }
    
    //boundary terms
    r0[nxSteps] = 0.5*(c/dx)*(u0[0] - u0[nxSteps-1]);
    r0[0] = 0.5*(c/dx)*(u0[1] - u0[nxSteps]);

    s0 = vd(nxSteps+1,0);
}

/*==========================================================================*/

void firstStep_leapfrog(vd &u0, vd &u1, vd &r0, vd &r1, vd &s0, vd &s1)
{
    //the 0.5 is because this first step is non-centred
    for(int i=1; i <= nxSteps; i++)
    {
        s1[i] = s0[i] + 0.5*alpha*( r0[i] - r0[i-1] );
    }
    //boundary term
    s1[0] = s0[0] + 0.5*alpha*( r0[0] - r0[nxSteps] );

    for(int i=0; i < nxSteps; i++)
    {
        r1[i] = r0[i] + alpha*( s1[i+1] - s1[i] );
    }
    //boundary term
    r1[nxSteps] = r0[nxSteps] + alpha*( s1[0] - s1[nxSteps] );

    for(int i=0; i <= nxSteps; i++)
    {
        u1[i] = u0[i] + dt*s1[i];
    }

    //rotate vectors 
    //the 0 quantities will be what we want i.e. values at current time
    u0.swap(u1);
    r0.swap(r1);
    s0.swap(s1);
}

/*==========================================================================*/

void oneStep_leapfrog(vd &u0, vd &u1, vd &r0, vd &r1, vd &s0, vd &s1)
{
    for(int i=1; i <= nxSteps; i++)
    {
        s1[i] = s0[i] + alpha*( r0[i] - r0[i-1] );
    }
    //boundary term    
    s1[0] = s0[0] + alpha*( r0[0] - r0[nxSteps] );

    for(int i=0; i < nxSteps; i++)
    {
        r1[i] = r0[i] + alpha*( s1[i+1] - s1[i] );
    }
    //boundary term
    r1[nxSteps] = r0[nxSteps] + alpha*( s1[0] - s1[nxSteps] );

    for(int i=0; i <= nxSteps; i++)
    {
        u1[i] = u0[i] + dt*s1[i];
    }

    //rotate vectors 
    //the 0 quantities will be what we want i.e. values at current time
    u0.swap(u1);
    r0.swap(r1);
    s0.swap(s1);
}

/*==========================================================================*/

void oneIter_icn(vd &r0, vd &r01, vd &r1, vd &s0, vd &s01, vd &s1)
{
    for(int i=1; i < nxSteps; i++)
    {
        r1[i] = r0[i] + 0.5*alpha*( s01[i+1] - s01[i-1] );
    }
    //boundary terms
    r1[nxSteps] = r0[nxSteps] + 0.5*alpha*( s01[0] - s01[nxSteps-1] );
    r1[0] = r0[0] + 0.5*alpha*( s01[1] - s01[nxSteps] );

    for(int i=1; i < nxSteps; i++)
    {
        s1[i] = s0[i] + 0.5*alpha*( r01[i+1] - r01[i-1] );
    }
    //boundary terms    
    s1[nxSteps] = s0[nxSteps] + 0.5*alpha*( r01[0] - r01[nxSteps-1] );
    s1[0] = s0[0] + 0.5*alpha*( r01[1] - r01[nxSteps] );
}

/*==========================================================================*/

void averageIter_icn(vd &r0, vd &r1, vd &s0, vd &s1)
{
    for(int i = 0; i <= nxSteps; i++)
    {
        r1[i] = 0.5*( r0[i] + r1[i] );
        s1[i] = 0.5*( s0[i] + s1[i] );
    }
}

/*==========================================================================*/

void calc_u_icn(vd &u0, vd &u1, vd &s0, vd &s1, vd &r0, vd &r1)
{
    for(int i=0; i<= nxSteps; i++)
    {
        u1[i] = u0[i] + 0.5*dt*( s0[i] + s1[i] );
    }
    
    //rotate vectors 
    //the 0 quantities are what we want i.e. values at current time
    u0.swap(u1);
    s0.swap(s1);
    r0.swap(r1);
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

    //this writes a script to animate in gnuplot
	animf << "set xrange [" << static_cast<int>(std::floor(xMin)) << ":" 
		<< static_cast<int>(std::ceil(xMax)) << "]\n";
	animf << "set yrange[" << -0.5*A << ":" << A+0.1 << "]\n";
	animf << "do for [i=0:" << ntSteps << "] {\n\ttime = " << dt << "*i\n";
	animf << "\ttitlevar = sprintf(\"time = %f\", time)\n";
	animf << "\tp 'data.dat' every :::i::i w l title titlevar\n";
	animf << "\tpause " << pause << "\n}\n";    
}