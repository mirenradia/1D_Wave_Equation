#include <iostream>
#include <fstream>
#include <ios>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib> //for exit()

typedef std::vector<double> vd;

struct GaussPars
{
    double a;
    double A;
    double x_0;
};

struct Pars
{
    double c;
    double dx;
    double dt;
    double alpha; //Courant factor
    double xMax;
    double xMin;
    double tMax;
    double pause;
    int nxSteps;
    int ntSteps;
    int n_icn;
    std::string alg;
};

//double c, dx, dt, alpha, xMax, xMin, tMax, a, A, x_0, pause;
//int nxSteps, ntSteps, n_icn;
//std::string alg;

void readPars(GaussPars&, Pars&, char*);
void init(GaussPars&, Pars&, vd&, vd&, vd&, vd&);
void init_icn(GaussPars&, Pars&, vd&, vd&, vd&, vd&);
void firstStep_leapfrog(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void oneStep_leapfrog(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void oneIter_icn(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void averageIter_icn(Pars&, vd&, vd&, vd&, vd&);
void calc_u_icn(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void writeData(Pars&, vd&, vd&, double);
void writeAnimationScript(GaussPars&, Pars&);

/*==========================================================================*/

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <parameter file>" << std::endl;
        exit(1);
    }

    GaussPars gp;
    Pars p;

    readPars(gp, p, argv[1]);

    //instantiate our vectors to hold solution and first derivatives
    //the 0 quantities are the values at the current time step
    //the 1 quantities are the values at the next time step
    vd x(p.nxSteps+1,0);
    vd u0(p.nxSteps+1,0);
    vd u1(p.nxSteps+1,0);
    vd r0(p.nxSteps+1,0);
    vd r1(p.nxSteps+1,0);
    vd s0(p.nxSteps+1,0);
    vd s1(p.nxSteps+1,0);
    double time {0.0};

    
    if(p.alg == "leapfrog")
    { 
        init(gp, p, u0, r0, s0, x);
        writeData(p, x, u0, time);
        firstStep_leapfrog(p, u0, u1, r0, r1, s0, s1);
        time += p.dt;
        writeData(p, x, u0, time);

        for(int t=1; t <= p.ntSteps; t++)
        {
            oneStep_leapfrog(p, u0, u1, r0, r1, s0, s1);
            time += p.dt;
            writeData(p, x, u0, time);
        }
	}
    else if(p.alg == "icn")
    {
        init_icn(gp, p, u0, r0, s0, x);
        writeData(p, x, u0, time);
        
        //instantiate 2 more vectors to hold intermediate values in algorithm
        //these will hold the "1/2" values in the ICN iterations
        vd r01(p.nxSteps+1,0);
        vd s01(p.nxSteps+1,0);

        for(int t=0; t <= p.ntSteps; t++)
        {
            for(int i_icn=0; i_icn <= p.n_icn; i_icn++)
            {
                if(i_icn == 0)
                {
                    oneIter_icn(p, r0, r0, r01, s0, s0, s01);
                }
                else if(i_icn != p.n_icn)
                {
                    averageIter_icn(p, r0, r01, s0, s01);
                    oneIter_icn(p, r0, r01, r1, s0, s01, s1);
                    
                    //swap vectors as this is not the final iteration
                    r1.swap(r01);
                    s1.swap(s01);
                }
                else
                {
                    averageIter_icn(p, r0, r01, s0, s01);
                    oneIter_icn(p, r0, r01, r1, s0, s01, s1);
                }
            }
            calc_u_icn(p, u0, u1, s0, s1, r0, r1);
            time += p.dt;
            writeData(p, x, u0, time);
        }
    }
    else
    {
        std::cerr << "In " << argv[1] << " there must be a line with either:\n";
        std::cerr << "alg = icn\nalg = leapfrog" << std::endl;
        return 1;
    }
    
	writeAnimationScript(gp, p);
    return 0;
}

/*==========================================================================*/

void readPars(GaussPars &gp, Pars &p, char* parf)
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
            if(line.substr(0,3) == "c =") {p.c = std::stod(line.substr(4));}
            if(line.substr(0,4) == "dt =") {p.dt = std::stod(line.substr(5));}
            if(line.substr(0,4) == "dx =") {p.dx = std::stod(line.substr(5));}
            if(line.substr(0,6) == "xMax =") {p.xMax = std::stod(line.substr(7));}
            if(line.substr(0,6) == "xMin =") {p.xMin = std::stod(line.substr(7));}
            if(line.substr(0,6) == "tMax =") {p.tMax = std::stod(line.substr(7));}
            if(line.substr(0,3) == "a =") {gp.a = std::stod(line.substr(4));}
            if(line.substr(0,5) == "x_0 =") {gp.x_0 = std::stod(line.substr(6));}
            if(line.substr(0,3) == "A =") {gp.A = std::stod(line.substr(4));}
            if(line.substr(0,5) == "alg =") {p.alg = line.substr(6);}
            if(line.substr(0,7) == "n_icn =") {p.n_icn = std::stoi(line.substr(8));}
            if(line.substr(0,7) == "pause =") {
            	p.pause = std::stod(line.substr(8));
            }

        }

        p.alpha = p.c*p.dt/p.dx;
        p.nxSteps = static_cast<int>((p.xMax-p.xMin)/p.dx);
        p.ntSteps = static_cast<int>(p.tMax/p.dt);
    }

/*==========================================================================*/

void init(GaussPars &gp, Pars &p, vd &u0, vd &r0, vd &s0, vd &x)
{
    for(int i=0; i <= p.nxSteps; i++)
    {
        x[i] = p.xMin + i*p.dx;
        u0[i] = gp.A*std::exp(-(x[i] - gp.x_0)*(x[i] - gp.x_0)/(2*gp.a*gp.a));
    }

    for(int i=0; i < p.nxSteps; i++)
    {
        r0[i] = (p.c/p.dx)*(u0[i+1] - u0[i]);
    }
    
    //boundary term
    r0[p.nxSteps] = (p.c/p.dx)*(u0[0] - u0[p.nxSteps]);

    s0 = vd(p.nxSteps+1,0);
}

/*==========================================================================*/

void init_icn(GaussPars &gp, Pars &p, vd &u0, vd &r0, vd &s0, vd &x)
{
    for(int i=0; i <= p.nxSteps; i++)
    {
        x[i] = p.xMin + i*p.dx;
        u0[i] = gp.A*std::exp(-(x[i] - gp.x_0)*(x[i] - gp.x_0)/(2*gp.a*gp.a));
    }

    for(int i=1; i < p.nxSteps; i++)
    {
        r0[i] = 0.5*(p.c/p.dx)*(u0[i+1] - u0[i-1]);
    }
    
    //boundary terms
    r0[p.nxSteps] = 0.5*(p.c/p.dx)*(u0[0] - u0[p.nxSteps-1]);
    r0[0] = 0.5*(p.c/p.dx)*(u0[1] - u0[p.nxSteps]);

    s0 = vd(p.nxSteps+1,0);
}

/*==========================================================================*/

void firstStep_leapfrog(Pars &p, vd &u0, vd &u1, vd &r0, vd &r1, vd &s0, vd &s1)
{
    //the 0.5 is because this first step is non-centred
    for(int i=1; i <= p.nxSteps; i++)
    {
        s1[i] = s0[i] + 0.5*p.alpha*( r0[i] - r0[i-1] );
    }
    //boundary term
    s1[0] = s0[0] + 0.5*p.alpha*( r0[0] - r0[p.nxSteps] );

    for(int i=0; i < p.nxSteps; i++)
    {
        r1[i] = r0[i] + p.alpha*( s1[i+1] - s1[i] );
    }
    //boundary term
    r1[p.nxSteps] = r0[p.nxSteps] + p.alpha*( s1[0] - s1[p.nxSteps] );

    for(int i=0; i <= p.nxSteps; i++)
    {
        u1[i] = u0[i] + p.dt*s1[i];
    }

    //swap vectors 
    //the 0 quantities will be what we want i.e. values at current time
    u0.swap(u1);
    r0.swap(r1);
    s0.swap(s1);
}

/*==========================================================================*/

void oneStep_leapfrog(Pars &p, vd &u0, vd &u1, vd &r0, vd &r1, vd &s0, vd &s1)
{
    for(int i=1; i <= p.nxSteps; i++)
    {
        s1[i] = s0[i] + p.alpha*( r0[i] - r0[i-1] );
    }
    //boundary term
    s1[0] = s0[0] + p.alpha*( r0[0] - r0[p.nxSteps] );

    for(int i=0; i < p.nxSteps; i++)
    {
        r1[i] = r0[i] + p.alpha*( s1[i+1] - s1[i] );
    }
    //boundary term
    r1[p.nxSteps] = r0[p.nxSteps] + p.alpha*( s1[0] - s1[p.nxSteps] );

    for(int i=0; i <= p.nxSteps; i++)
    {
        u1[i] = u0[i] + p.dt*s1[i];
    }

    //swap vectors 
    //the 0 quantities will be what we want i.e. values at current time
    u0.swap(u1);
    r0.swap(r1);
    s0.swap(s1);
}

/*==========================================================================*/

void oneIter_icn(Pars &p, vd &r0, vd &r01, vd &r1, vd &s0, vd &s01, vd &s1)
{
    for(int i=1; i < p.nxSteps; i++)
    {
        r1[i] = r0[i] + 0.5*p.alpha*( s01[i+1] - s01[i-1] );
    }
    //boundary terms
    r1[p.nxSteps] = r0[p.nxSteps] + 0.5*p.alpha*( s01[0] - s01[p.nxSteps-1] );
    r1[0] = r0[0] + 0.5*p.alpha*( s01[1] - s01[p.nxSteps] );

    for(int i=1; i < p.nxSteps; i++)
    {
        s1[i] = s0[i] + 0.5*p.alpha*( r01[i+1] - r01[i-1] );
    }
    //boundary terms    
    s1[p.nxSteps] = s0[p.nxSteps] + 0.5*p.alpha*( r01[0] - r01[p.nxSteps-1] );
    s1[0] = s0[0] + 0.5*p.alpha*( r01[1] - r01[p.nxSteps] );
}

/*==========================================================================*/

void averageIter_icn(Pars &p, vd &r0, vd &r1, vd &s0, vd &s1)
{
    for(int i = 0; i <= p.nxSteps; i++)
    {
        r1[i] = 0.5*( r0[i] + r1[i] );
        s1[i] = 0.5*( s0[i] + s1[i] );
    }
}

/*==========================================================================*/

void calc_u_icn(Pars &p, vd &u0, vd &u1, vd &s0, vd &s1, vd &r0, vd &r1)
{
    for(int i=0; i<= p.nxSteps; i++)
    {
        u1[i] = u0[i] + 0.5*p.dt*( s0[i] + s1[i] );
    }
    
    //swap vectors 
    //the 0 quantities are what we want i.e. values at current time
    u0.swap(u1);
    s0.swap(s1);
    r0.swap(r1);
}

/*==========================================================================*/

void writeData(Pars &p, vd &x, vd &u0, double time)
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

	for(int i = 0; i <= p.nxSteps; i++)
	{
		dataf << x[i] << "\t" << u0[i] << "\n";
	}

	dataf << "\n";
}

/*==========================================================================*/

void writeAnimationScript(GaussPars &gp, Pars &p)
{
	std::ofstream animf("animation.gpi");

	if(!animf)
	{
		std::cerr << "Uh oh, could not open animation.gpi for writing" 
				<< std::endl;
		return;
	}

    //this writes a script to animate in gnuplot
	animf << "set xrange [" << static_cast<int>(std::floor(p.xMin)) << ":" 
		<< static_cast<int>(std::ceil(p.xMax)) << "]\n";
	animf << "set yrange[" << -0.5*gp.A << ":" << gp.A+0.1 << "]\n";
	animf << "do for [i=0:" << p.ntSteps << "] {\n\ttime = " << p.dt << "*i\n";
	animf << "\ttitlevar = sprintf(\"time = %f\", time)\n";
	animf << "\tp 'data.dat' every :::i::i w l title titlevar\n";
	animf << "\tpause " << p.pause << "\n}\n";    
}