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
    int tSkip;
    int xSkip;
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

void readPars(GaussPars&, Pars&, char*);
void init_leapfrog(GaussPars&, Pars&, vd&, vd&, vd&, vd&);
void init_icn(GaussPars&, Pars&, vd&, vd&, vd&, vd&);
void firstStep_leapfrog(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void oneStep_leapfrog(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void oneIter_icn(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void averageIter_icn(Pars&, vd&, vd&, vd&, vd&);
void calc_u_icn(Pars&, vd&, vd&, vd&, vd&, vd&, vd&);
void oneStep_icn(Pars&, vd&, vd&, vd&, vd&, vd&, vd&, vd&, vd&);
void writeData(Pars&, vd&, vd&, double, std::string);
void writeAnimationScript(GaussPars&, Pars&, std::string, std::string);
void leapfrogMainLoop(GaussPars&, Pars&);
void icnMainLoop(GaussPars&, Pars&);
void leapfrogMainLoop_conv_test(GaussPars&, Pars&);
void icnMainLoop_conv_test(GaussPars&, Pars&);


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
    
    if(p.alg == "leapfrog")
    { 
        leapfrogMainLoop(gp, p);
	}
    else if(p.alg == "icn")
    {
        icnMainLoop(gp, p);
    }
    else if(p.alg == "leapfrog conv test")
    {
        leapfrogMainLoop_conv_test(gp, p);
    }
    else if(p.alg == "icn conv test")
    {
    	icnMainLoop_conv_test(gp, p);
    }
    else
    {
        std::cerr << "In " << argv[1] << " there must be a line with either:\n";
        std::cerr << "alg = icn\n";
        std::cerr << "alg = leapfrog\n";
        std::cerr << "alg = leapfrog conv test\n";
        std::cerr << "alg = icn conv test" << std::endl;
        return 1;
    }
    
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
        if(line.substr(0,7) == "tSkip =") {p.tSkip = std::stoi(line.substr(8));}
        if(line.substr(0,7) == "xSkip =") {p.xSkip = std::stoi(line.substr(8));}
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

void printPars(GaussPars &gp, Pars &p)
{
	std::cout << std::setprecision(15);

	std::cout << "Algorithm Parameters:\n";
	std::cout << "c       = " << p.c << "\n";
	std::cout << "dt      = " << p.dt << "\n";
	std::cout << "dx      = " << p.dx << "\n";
    std::cout << "tSkip   = " << p.tSkip << "\n";
    std::cout << "xSkip   = " << p.xSkip << "\n";
	std::cout << "alpha   = " << p.alpha << "\n";
	std::cout << "xMin    = " << p.xMin << "\n";
	std::cout << "xMax    = " << p.xMax << "\n";
	std::cout << "tMax    = " << p.tMax << "\n";
	std::cout << "nxSteps = " << p.nxSteps << "\n";
	std::cout << "ntSteps = " << p.ntSteps << "\n";
	std::cout << "alg     = " << p.alg << "\n";
	if(!(p.alg.find("icn") == std::string::npos))
		std::cout << "n_icn   = " << p.n_icn << "\n";

	std::cout << "\nInitial Gaussian Parameters:\n";
	std::cout << "a       = " << gp.a << "\n";
	std::cout << "x_0     = " << gp.x_0 << "\n";
	std::cout << "A       = " << gp.A << "\n";
}

/*==========================================================================*/

void init_leapfrog(GaussPars &gp, Pars &p, vd &u0, vd &r0, vd &s0, vd &x)
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

    for(int i=0; i <= p.nxSteps; i++)
    {
        r0[i] = -(x[i] - gp.x_0)*u0[i]/(gp.a*gp.a);
    }
    
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

void oneStep_icn(Pars &p, vd &u0, vd &u1, vd &r0, vd &r01, vd &r1, vd &s0, 
	 vd &s01, vd &s1)
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
}

/*==========================================================================*/

void writeData(Pars &p, vd &x, vd &u0, double time, std::string filename)
{
	std::ofstream dataf(filename, (time == 0) ? std::ios::out : std::ios::app);

	if(!dataf)
	{
		std::cerr << "Uh oh, could not open" << filename << "for writing" << 
        std::endl;
		exit(1);
	}
	dataf << std::scientific << std::setprecision(14);
	dataf << "#time = " << time << "\n";

	for(int i = 0; i <= p.nxSteps/p.xSkip; i++)
	{
		dataf << x[p.xSkip*i] << "\t" << u0[p.xSkip*i] << "\n";
	}

	dataf << "\n";
}

/*==========================================================================*/

void writeData(Pars &p, vd &x, vd &u0, vd &u1, double time, 
    std::string filename)
{
    std::ofstream dataf(filename, (time == 0) ? std::ios::out : std::ios::app);

    if(!dataf)
    {
        std::cerr << "Uh oh, could not open" << filename << "for writing" << 
        std::endl;
        exit(1);
    }
    dataf << std::scientific << std::setprecision(14);
    dataf << "#time = " << time << "\n";

    for(int i = 0; i <= p.nxSteps/p.xSkip; i++)
    {
        dataf << x[p.xSkip*i] << "\t" << u0[p.xSkip*i] << "\t" 
        << u1[p.xSkip*i] << "\n";
    }

    dataf << "\n";
}

/*==========================================================================*/

void writeAnimationScript(GaussPars &gp, Pars &p, std::string datafname, 
    std::string animfname)
{
	std::ofstream animf(animfname);

	if(!animf)
	{
		std::cerr << "Uh oh, could not open animation.gpi for writing" 
				<< std::endl;
		return;
	}

    //this writes a script to animate in gnuplot
	animf << "set xrange [" << static_cast<int>(std::floor(p.xMin)) << ":" 
		<< static_cast<int>(std::ceil(p.xMax)) << "]\n";
	if(p.alg.find("conv test") == std::string::npos)
    {
        animf << "set yrange[" << -0.5*gp.A << ":" << gp.A+0.1 << "]\n";
	}
    animf << "do for [i=0:" << p.ntSteps << "] {\n\ttime = " << p.dt << "*i\n";
	if(p.alg.find("conv test") == std::string::npos)
    {
        animf << "\ttitlevar = sprintf(\"time = %f\", time)\n";
        animf << "\tp '" << datafname << "' every :::i::i w l title titlevar\n";
        animf << "\tpause " << p.pause << "\n}\n";    
    }
    else
    {
        animf << "\tset title sprintf(\"time = %f\", time)\n";
        animf << "\tp '" << datafname << "' using 1:2 every :::i::i w l";
        animf << " title '|u_c-u_m|', \\\n";
        animf << "\t'" << datafname << "' using 1:3 every :::i::i w l";
        animf << " title '|u_m-u_f|',\n";
        animf << "\tpause " << p.pause << "\n}\n";
    }
}

/*==========================================================================*/

void leapfrogMainLoop(GaussPars &gp, Pars &p)
{
    printPars(gp, p);

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

    init_leapfrog(gp, p, u0, r0, s0, x);
    writeData(p, x, u0, time, "data.dat");
    firstStep_leapfrog(p, u0, u1, r0, r1, s0, s1);
    time += p.dt;
    
    if(p.tSkip == 1)
    {
        writeData(p, x, u0, time, "data.dat");
    }

    for(int t=2; t <= p.ntSteps; t++)
    {
        oneStep_leapfrog(p, u0, u1, r0, r1, s0, s1);
        time += p.dt;
        if(t % p.tSkip == 0)
        {
            writeData(p, x, u0, time, "data.dat");
        }
    }

    writeAnimationScript(gp, p, "data.dat", "animation.gpi");
}

/*==========================================================================*/

void icnMainLoop(GaussPars &gp, Pars &p)
{
    printPars(gp, p);

    //instantiate our vectors to hold solution and first derivatives
    //the 0 quantities are the values at the current time step
    //the 1 quantities are the values at the next time step
    //the 01 quantities will hold the intermediate "1/2" values in the
    //ICN algorithm
    vd x(p.nxSteps+1,0);
    vd u0(p.nxSteps+1,0);
    vd u1(p.nxSteps+1,0);
    vd r0(p.nxSteps+1,0);
    vd r01(p.nxSteps+1,0);
    vd r1(p.nxSteps+1,0);
    vd s0(p.nxSteps+1,0);
    vd s01(p.nxSteps+1,0);
    vd s1(p.nxSteps+1,0);
    double time {0.0};
    
    init_icn(gp, p, u0, r0, s0, x);
    writeData(p, x, u0, time, "data.dat");

    for(int t=1; t <= p.ntSteps; t++)
    {
    	oneStep_icn(p, u0, u1, r0, r01, r1, s0, s01, s1);
        time += p.dt;
        if(t % p.tSkip == 0)
        {
            writeData(p, x, u0, time, "data.dat");
        }
    }   

    writeAnimationScript(gp, p, "data.dat", "animation.gpi");
}

/*==========================================================================*/

void leapfrogMainLoop_conv_test(GaussPars &gp, Pars &p_c)
{
    //the input parameters will be for the coarse grid
    //create 2 new Parameter structs for the medium and fine grids
    //In this function:
    //_c denote coarse grid values
    //_m denote medium grid values
    //_f denote fine grid values
    Pars p_m {p_c};
    Pars p_f {p_c};

    //refine grid size by 2 for medium and fine grids
    p_m.dt *= 0.5;
    p_m.dx *= 0.5;
    p_m.ntSteps *= 2;
    p_m.nxSteps *= 2;
    p_f.dt *= 0.25;
    p_f.dx *= 0.25;
    p_f.ntSteps *= 4;
    p_f.nxSteps *= 4;

    std::cout << "\nCoarse ";
    printPars(gp, p_c);
    std::cout << "\nMedium ";
    printPars(gp, p_m);
    std::cout << "\nFine ";
    printPars(gp, p_f);


    //instantiate our vectors to hold solution and first derivatives
    //the 0 quantities are the values at the current time step
    //the 1 quantities are the values at the next time step
    vd x_c(p_c.nxSteps+1,0);
    vd u0_c(p_c.nxSteps+1,0);
    vd u1_c(p_c.nxSteps+1,0);
    vd r0_c(p_c.nxSteps+1,0);
    vd r1_c(p_c.nxSteps+1,0);
    vd s0_c(p_c.nxSteps+1,0);
    vd s1_c(p_c.nxSteps+1,0);

    vd x_m(p_m.nxSteps+1,0);
    vd u0_m(p_m.nxSteps+1,0);
    vd u1_m(p_m.nxSteps+1,0);
    vd r0_m(p_m.nxSteps+1,0);
    vd r1_m(p_m.nxSteps+1,0);
    vd s0_m(p_m.nxSteps+1,0);
    vd s1_m(p_m.nxSteps+1,0);

    vd x_f(p_f.nxSteps+1,0);
    vd u0_f(p_f.nxSteps+1,0);
    vd u1_f(p_f.nxSteps+1,0);
    vd r0_f(p_f.nxSteps+1,0);
    vd r1_f(p_f.nxSteps+1,0);
    vd s0_f(p_f.nxSteps+1,0);
    vd s1_f(p_f.nxSteps+1,0);

    double time {0.0};

    //instantiate two vectors to hold the difference between solutions
    vd diff_c_m(p_c.nxSteps+1,0);
    vd diff_m_f(p_c.nxSteps+1,0);

    init_leapfrog(gp, p_c, u0_c, r0_c, s0_c, x_c);
    init_leapfrog(gp, p_m, u0_m, r0_m, s0_m, x_m);
    init_leapfrog(gp, p_f, u0_f, r0_f, s0_f, x_f);

    for(int i=0; i <= p_c.nxSteps; i++)
    {
        diff_c_m[i] = std::abs(u0_c[i] - u0_m[2*i]);
        diff_m_f[i] = std::abs(u0_m[2*i] - u0_f[4*i]);
    }

    writeData(p_c, x_c, diff_c_m, diff_m_f, time, "conv_test_data.dat");

    firstStep_leapfrog(p_c, u0_c, u1_c, r0_c, r1_c, s0_c, s1_c);
    firstStep_leapfrog(p_m, u0_m, u1_m, r0_m, r1_m, s0_m, s1_m);
    firstStep_leapfrog(p_f, u0_f, u1_f, r0_f, r1_f, s0_f, s1_f);

    for(int t=2; t <= p_f.ntSteps; t++)
    {
        time += p_f.dt;
        oneStep_leapfrog(p_f, u0_f, u1_f, r0_f, r1_f, s0_f, s1_f);
        
        if(t % 2 == 0)
        {
            oneStep_leapfrog(p_m, u0_m, u1_m, r0_m, r1_m, s0_m, s1_m);
        }    
        
        if(t % 4 == 0)
        {
            oneStep_leapfrog(p_c, u0_c, u1_c, r0_c, r1_c, s0_c, s1_c);

            for(int i=0; i <= p_c.nxSteps; i++)
            {
                diff_c_m[i] = std::abs(u0_c[i] - u0_m[2*i]);
                diff_m_f[i] = std::abs(u0_m[2*i] - u0_f[4*i]);
            }
            
            if(t % (4*p_c.tSkip) == 0)
            {
                writeData(p_c, x_c, diff_c_m, diff_m_f, time, 
                    "conv_test_data.dat");
            }
        }
    }
    writeAnimationScript(gp, p_c, "conv_test_data.dat", 
        "conv_test_animation.gpi");
}

/*==========================================================================*/

void icnMainLoop_conv_test(GaussPars &gp, Pars &p_c)
{
    //the input parameters will be for the coarse grid
    //create 2 new Parameter structs for the medium and fine grids
    //In this function:
    //_c denote coarse grid values
    //_m denote medium grid values
    //_f denote fine grid values
    Pars p_m {p_c};
    Pars p_f {p_c};

    //refine grid size by 2 for medium and fine grids
    p_m.dt *= 0.5;
    p_m.dx *= 0.5;
    p_m.ntSteps *= 2;
    p_m.nxSteps *= 2;
    p_f.dt *= 0.25;
    p_f.dx *= 0.25;
    p_f.ntSteps *= 4;
    p_f.nxSteps *= 4;

    std::cout << "\nCoarse ";
    printPars(gp, p_c);
    std::cout << "\nMedium ";
    printPars(gp, p_m);
    std::cout << "\nFine ";
    printPars(gp, p_f);


    //instantiate our vectors to hold solution and first derivatives
    //the 0 quantities are the values at the current time step
    //the 1 quantities are the values at the next time step
    //the 01 quantities will hold the intermediate "1/2" values in the
    //ICN algorithm
    vd x_c(p_c.nxSteps+1,0);
    vd u0_c(p_c.nxSteps+1,0);
    vd u1_c(p_c.nxSteps+1,0);
    vd r0_c(p_c.nxSteps+1,0);
    vd r01_c(p_c.nxSteps+1,0);
    vd r1_c(p_c.nxSteps+1,0);
    vd s0_c(p_c.nxSteps+1,0);
    vd s01_c(p_c.nxSteps+1,0);
    vd s1_c(p_c.nxSteps+1,0);

    vd x_m(p_m.nxSteps+1,0);
    vd u0_m(p_m.nxSteps+1,0);
    vd u1_m(p_m.nxSteps+1,0);
    vd r0_m(p_m.nxSteps+1,0);
    vd r01_m(p_m.nxSteps+1,0);
    vd r1_m(p_m.nxSteps+1,0);
    vd s0_m(p_m.nxSteps+1,0);
    vd s01_m(p_m.nxSteps+1,0);
    vd s1_m(p_m.nxSteps+1,0);

    vd x_f(p_f.nxSteps+1,0);
    vd u0_f(p_f.nxSteps+1,0);
    vd u1_f(p_f.nxSteps+1,0);
    vd r0_f(p_f.nxSteps+1,0);
    vd r01_f(p_f.nxSteps+1,0);
    vd r1_f(p_f.nxSteps+1,0);
    vd s0_f(p_f.nxSteps+1,0);
    vd s01_f(p_f.nxSteps+1,0);
    vd s1_f(p_f.nxSteps+1,0);

    double time {0.0};

    //instantiate two vectors to hold the difference between solutions
    vd diff_c_m(p_c.nxSteps+1,0);
    vd diff_m_f(p_c.nxSteps+1,0);

    init_icn(gp, p_c, u0_c, r0_c, s0_c, x_c);
    init_icn(gp, p_m, u0_m, r0_m, s0_m, x_m);
    init_icn(gp, p_f, u0_f, r0_f, s0_f, x_f);

    for(int i=0; i <= p_c.nxSteps; i++)
    {
        diff_c_m[i] = std::abs(u0_c[i] - u0_m[2*i]);
        diff_m_f[i] = std::abs(u0_m[2*i] - u0_f[4*i]);
    }

    writeData(p_c, x_c, diff_c_m, diff_m_f, time, "conv_test_data.dat");

    for(int t=1; t<= p_f.ntSteps; t++)
    {
        oneStep_icn(p_f, u0_f, u1_f, r0_f, r01_f, r1_f, s0_f, s01_f, s1_f);
        
        if(t % 2 == 0)
        {
            oneStep_icn(p_m, u0_m, u1_m, r0_m, r01_m, r1_m, s0_m, s01_m, s1_m);
        }    
        
        if(t % 4 == 0)
        {
            oneStep_icn(p_c, u0_c, u1_c, r0_c, r01_c, r1_c, s0_c, s01_c, s1_c);

            for(int i=0; i <= p_c.nxSteps; i++)
            {
                diff_c_m[i] = std::abs(u0_c[i] - u0_m[2*i]);
                diff_m_f[i] = std::abs(u0_m[2*i] - u0_f[4*i]);
            }

            if(t % (4*p_c.tSkip) == 0)
            {
                writeData(p_c, x_c, diff_c_m, diff_m_f, time, 
                    "conv_test_data.dat");
            }
        }
        time += p_f.dt;
    }

    writeAnimationScript(gp, p_c, "conv_test_data.dat", 
        "conv_test_animation.gpi");
}
