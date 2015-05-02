// particles_v2.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "gnuplot_i.hpp"
#include <iostream>
#include <memory>

#define _USE_MATH_DEFINES

#include "math.h"
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include <fstream>
#include <ctime> 


// m =91aem q =1 2 3


//double m_PI = 3.1415;
double m = 2;// 91 * 1.66e-27;
double m_e = 2;// 9.109e-28;// mass
double dpi = 2*M_PI;
double k = 1;// 1.380648813e-23;
double T = 1;// 3600;
double q1 = 1;// (1.6e-23);
double q2 = 2;// *(1.6e-23);
double q3 = 4;// *(1.6e-23);
double kc = 1;//(1./(4*M_PI*8.8554187817e-12));


/***********************************/
int CP = 10;//count of particles		
int N = 1000;//time steps
double critical_rad = 110.8;
/***********************************/

class vec3d
{
public:

	double x, y, z;

	vec3d(){
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	vec3d(double _x, double _y, double _z){
		this->x = _x;
		this->y = _y;
		this->z = _z;
	}

};


vec3d operator*(double a, vec3d b)
{
	return vec3d(a*b.x, a*b.y, a*b.z);
}

vec3d operator%(vec3d a, vec3d b)
{
	return vec3d(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

vec3d operator+(vec3d a, vec3d b)
{
	return vec3d(a.x + b.x, a.y + b.y, a.z + b.z);
}

vec3d operator-(vec3d a, vec3d b)
{
	return vec3d(a.x - b.x, a.y - b.y, a.z - b.z);
}

double operator*(vec3d a, vec3d b)
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}




template <typename T>
T random(T Min, T Max) {
	return Min + static_cast <T> (rand()) / (static_cast <T> (RAND_MAX / (Max - Min)));
}





double Gauss(vec3d v, double sigma)
{
	return pow((2 * 1 / (sigma*sqrt(2 * M_PI))), 3)*exp(-(v*v) / (2 * sigma*sigma));
}

double Maxwell(vec3d vr, double sigma)
{
	return sqrt(2 * m / dpi*sigma)*exp(-(m*vr*vr) / (2 * k*T));
}





vec3d GenerateRM(double sigma) {
	double dr, dA, dB, dp;
	do {
		//dr = random(0.0,1.0*sigma);
		double U = random(0.0, 1.0);
		double V = random(0.0, 1.0);
		dA = M_PI*U;
		dB = acos(2 * V - 1);
		dp = random(0.0, sigma);
	} while (dp < Maxwell(vec3d(sin(dB)*cos(dA), sin(dB)*sin(dA), cos(dB)), sigma));
	return dp*vec3d(cos(dB), sin(dB)*cos(dA), sin(dB)*sin(dA));//vec3d(sin(dB)*cos(dA),sin(dB)*sin(dA),cos(dB));
}

vec3d GenerateR(double sigma, double length) {
	double dr, dA, dB;
	double A = random(0.0, M_PI*2.0);
	double L = random(0.0, length);
	double G = random(0.0, 1.0);
	dr = sqrtl(-2 * logl(G))*sigma;
	return vec3d(dr*cos(A), dr*sin(A), L);
}


double module_btw2(vec3d first, vec3d second)
{
	return sqrt(pow((first.x - second.x), 2.) + pow((first.y - second.y), 2.) + pow((first.z - second.z), 2.));
}

vec3d pwr_btw2(vec3d first, vec3d second, double QQ)
{
	vec3d r = first - second;
	double yy2 = pow(module_btw2(first, second), 3);
	double yy = pow(QQ, 2);
	double yy1 = yy*kc;
	double yy3 = (yy1 / yy2);

	if (yy2 == 0){
		return 0 * r;
	}
	return yy3*r;
}

double absolute(vec3d G)
{
	return sqrt(G.x*G.x + G.y*G.y + G.z*G.z);
}

bool _contains(std::vector<int> tar, int cp)
{
	for (int ss = 0; ss< tar.size(); ss++)
	{
		if (tar[ss] == cp){
			return true;
		}
	}
	return false;
}

vec3d CuloPower(vec3d p1, vec3d p2, double QQ)
{
	double module = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	//double QQ = p1.q*p2.q;
	vec3d temp = p2 - p1;
	vec3d yy = (1 / pow(module, 3))*temp;
	if (module <= 1e-16)
	{
		return vec3d(0, 0, 0);
	}
	return (QQ / kc*yy);
}


int main()
{
	std::ofstream FinRx;
	FinRx.open("inRx1.txt");
	std::ofstream FinVx;
	FinVx.open("inVx1.txt");
	std::ofstream FoutRx;
	FoutRx.open("outRx1.txt");
	std::ofstream FoutVx;
	FoutVx.open("outVx1.txt");

	double _q[CP];
	double U = 10000;
	int mn;
	double dT = 0.001;

	std::vector<vec3d> inR_vect(CP);
	std::vector<vec3d> outR_vect(CP);
	std::vector<vec3d> inV_vect(CP);
	std::vector<vec3d> outV_vect(CP);
	std::vector<vec3d> p_frc(CP);

	vec3d B(1, 0, 0);
	vec3d E(0, 0, 10000);
	double L = 1;//10 cm

	
	// generate beam
	for (int gg = 0; gg < CP; gg++)
	{
		double q = 1;
		_q[gg] = q;
		inR_vect[gg] = GenerateR(10.2, 10.5);//Height Lenght
		inV_vect[gg] = GenerateRM(k*T);
		if (gg>80)
			q = q2;
		if (gg > 120)
			q = q3;
	}

	//calc particles
	for (int i = 0; i < N; i++)
	{
		std::cout << i << '\n';
		vec3d temp(0, 0, 0);
		for (int gg2 = 0; gg2 < CP; gg2++) {

			for (int p = 0; p < CP; p++) {
				if (p != gg2) {
					temp = temp + CuloPower(outR_vect[gg2], outR_vect[p], _q[p] * _q[gg2]);
				}
			}
			//vec3d ez(0, 0, 1);
			if (outR_vect[gg2].z < L) // electric field
			{
				p_frc[gg2] = temp + (_q[gg2] / L) * E;
				outV_vect[gg2] = outV_vect[gg2] + dT*p_frc[gg2];
				outR_vect[gg2] = outR_vect[gg2] + (dT * outV_vect[gg2]);

				FinRx << outR_vect[gg2].x << " ";
				FinRx << outR_vect[gg2].y << " ";
				FinRx << outR_vect[gg2].z << "\n";

				FinVx << outV_vect[gg2].x << " ";
				FinVx << outV_vect[gg2].y << " ";
				FinVx << outV_vect[gg2].z << "\n";


			}

			// magnetic field
			if (outR_vect[gg2].z >= L)
			{
				p_frc[gg2] = temp + _q[gg2] * (outV_vect[gg2] % B);
				outV_vect[gg2] = outV_vect[gg2] + (1 / m_e) * (dT * p_frc[gg2]);
				outR_vect[gg2] = outR_vect[gg2] + (dT * outV_vect[gg2]);

				/*FinRx << inR_vect[gg2].x << " ";
				FinRx << inR_vect[gg2].y << " ";
				FinRx << inR_vect[gg2].z << "\n";*/

				//				FinRx << inV_vect[gg2].x / inV_vect[gg2].z << " ";
				//				FinRx << inR_vect[gg2].y << " ";
				//				FinRx << 0 << "\n";
				//
				//				FinVx << inV_vect[gg2].x << " ";
				//				FinVx << inV_vect[gg2].y << " ";
				//				FinVx << inV_vect[gg2].z << "\n";

				/*FinRx << outR_vect[gg2].x << " ";
				FinRx << outR_vect[gg2].y << " ";
				FinRx << outR_vect[gg2].z << "\n";

				FinVx << outV_vect[gg2].x << " ";
				FinVx << outV_vect[gg2].y << " ";
				FinVx << outV_vect[gg2].z << "\n";*/
			}

		}
	}

	FinRx.close();
	FinVx.close();
	FoutRx.close();
	FoutVx.close();
	Gnuplot gp1;
	//	gp1.cmd("splot 'inRx1.txt' with dots");
	//	getchar();
	gp1.cmd("splot 'inRx1.txt' with dots");
	getchar();
	return 0;
}


