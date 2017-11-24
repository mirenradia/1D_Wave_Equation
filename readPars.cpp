#include <ifstream>
#include <iostream>
#include <string>

void readPars(double &c, double &dt, double &dx, double &alpha, double &xMax,
             double &xMin, int &nxSteps, double &tMax, int &ntSteps,
             double &a, double &A, double &x_0, char* parf)
    {
        std::ifstream parstream(parf)

        if (!parstream)
        {
            std::cerr << "Uh oh, could not open " << parf << " for reading." << std::endl;
            exit(1)
        }

        std::string line

        while(parstream)
        {
            std::getline(parstream, line)
            if(line.substr(0,3) == "c =") {c = std::stod(line.substr(4));}
            if(line.substr(0,4) == "dt =") {dt = std::stod(line.substr(5));}
            if(line.substr(0,4) == "dx =") {dx = std::stod(line.substr(5));}
            if(line.substr(0,6) == "xMax =") {xMax = std::stod(line.substr(7));}
            if(line.substr(0,6) == "xMin =") {xMin = std::stod(line.substr(7));}
            if(line.substr(0,6) == "tMax =") {tMax = std::stod(line.substr(7));}
            if(line.substr(0,3) == "a =") {a = std::stod(line.substr(4));}
            if(line.substr(0,5) == "x_0 =") {x_0 = std::stod(line.substr(6));}
        }

        alpha = c*dt/dx;
        nxSteps = static_cast<int>((xMax-xMin)/dx);
        ntSteps = static_cast<int>(tMax/dt);
    }
