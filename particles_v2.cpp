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
double m = 2* 1.66e-27;
//double m_e =  9.109e-28;// mass
double dpi = 2*M_PI;
double k = 0.00001;// 1.380648813e-23;
double T =  3600;
double q1 = 1*(1.6e-19);
double q2 = 2 *(1.6e-19);
double q3 = 4 *(1.6e-19);
double kc = (1./(4*M_PI*(8.8554187817e-12)));

double tangens30 = 1 / sqrt(3);
/***********************************/
int CP = 500;//count of particles		
int N = 15000;//time steps
double L = 0.1;//10 cm
double R0 = 0.33; // 33 cm
double Width = 0.015;//1.5 cm
//double critical_rad = 110.8;
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
		//double n1=random(0.0, 1.0);
		//double n2=random(0.0, 1.0);
		//double n3=random(0.0, 1.0);

		return vec3d(0, 0, 0);
	}
	return (QQ / kc*yy);
}


double intersection(double inputZ)
{
	double temp = (-tangens30*(inputZ - L) + R0);
	return temp;
}
double intersection_INPUT(double inputZ)
{
	double temp = (-tangens30*(inputZ - L) + R0);
	return temp;
}
double PIFAGOR(double z,double y)
{
	double temp = sqrt(z*z+y*y) ;
	return temp;
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
	std::vector<vec3d> inR_vect(CP);
	std::vector<vec3d> outR_vect(CP);
	std::vector<vec3d> inV_vect(CP);
	std::vector<vec3d> outV_vect(CP);
	std::vector<vec3d> p_frc(CP);

	double _q[CP];
	double _m[CP];
	double U = 10000;
	int mn;
	double dT = 0.00000001;//10e-8

	       dT = 0.000000001;

	vec3d B(0.0811, 0, 0);
	vec3d E(0, 0, 100000);
	//double L = 0.1;//10 cm
	bool leftE[CP];
	
	// generate beam
	for (int gg = 0; gg < CP; gg++)
	{
		_q[gg] = q1;
		_m[gg] = m;
		inR_vect[gg] = GenerateR(0.0075, 0.03);//Height Lenght
		inV_vect[gg] = GenerateRM(1*k*T);
		if (gg>500)
			_q[gg] = q2;
		_m[gg] = 2*m;
		if (gg>8000)
			_q[gg] = q3;
		_m[gg] = m;
		if (gg>10000)
			_q[gg] = 2*q1;
		_m[gg] = 3*m;
		
		
		

		/*FinRx << inV_vect[gg].x << " ";
	    FinRx << inV_vect[gg].y << " ";
		FinRx << inV_vect[gg].z << "\n";

		FoutRx << inR_vect[gg].x << " ";
		FoutRx << inR_vect[gg].y << " ";
		FoutRx << inR_vect[gg].z << "\n";*/
	}
	
	for (int gg_gg = 0; gg_gg < CP; gg_gg++)
	{
		outV_vect[gg_gg] = inV_vect[gg_gg];
		outR_vect[gg_gg] = inR_vect[gg_gg];
		leftE[gg_gg] = false;
	}

	//calc particles
	
	for (int i = 0; i < N; i++)
	{
		//////////////
		//continue;/////TURN OFF CALC
		//////////////
		std::cout << i << '\n';
		for (int gg2 = 0; gg2 < CP; gg2++) 
		{
			//////// turn on if no need B
			if ((outR_vect[gg2].z) > L)
				continue;
			///////////
			//here we can exept particles///////////////////////////
			//if 
			bool ccc0 = (outR_vect[gg2].y > intersection(outR_vect[gg2].z));
			bool ccc1 = ((PIFAGOR(outR_vect[gg2].y - Width-R0, outR_vect[gg2].z - L) < (R0 - Width)));
			bool ccc2 = ((PIFAGOR(outR_vect[gg2].y + Width-R0, outR_vect[gg2].z - L) > (R0 + Width)));
			bool ccc =  (ccc1||ccc2);
			if (ccc0 || ccc)
			continue;	
			
			//
			////////////////////////////////////////////////////////
			vec3d temp(0, 0, 0);
			for (int p = 0; p < CP; p++) 
			{
				if (p != gg2) 
				{
					temp = temp + CuloPower(outR_vect[gg2], outR_vect[p], _q[p] * _q[gg2]);
				}
			}

			if (outR_vect[gg2].z <= L) // electric field
			{
				
				p_frc[gg2] = temp + ((_q[gg2]) * E) ;
				outV_vect[gg2] = outV_vect[gg2] + (1. / _m[gg2])*dT*p_frc[gg2];
				outR_vect[gg2] = outR_vect[gg2] + (dT * outV_vect[gg2]);// +(1. / m) *(dT*dT / 2)*(p_frc[gg2]);
				
				FinRx << outR_vect[gg2].x << " ";
				FinRx << outR_vect[gg2].y << " ";
			    FinRx << outR_vect[gg2].z << "\n";

				//FinVx << outV_vect[gg2].x << " ";
				//FinVx << outV_vect[gg2].y << " ";
				//FinVx << outV_vect[gg2].z << "\n";
			}
			
			

			// magnetic field
			if ((outR_vect[gg2].z) > L)
			{ 
				
				if (!leftE[gg2])
				{

					//FinVx << outR_vect[gg2].x / outR_vect[gg2].z << " ";
					//FinVx << outR_vect[gg2].x << " ";
					//FinVx << 0 << "\n";
					//leftE[gg2] = true;

					//FinRx << outR_vect[gg2].x << " ";
					//FinRx << outR_vect[gg2].y << " ";
					//FinRx << outR_vect[gg2].z << "\n";
				}


				p_frc[gg2] = temp + _q[gg2] * (outV_vect[gg2] % B);// -((_m[gg2] / R0)* outV_vect[gg2]);
				outV_vect[gg2] = outV_vect[gg2] + (1. / _m[gg2]) * (dT * p_frc[gg2]);
				outR_vect[gg2] = outR_vect[gg2] + (dT * outV_vect[gg2]);// +(1. / m) * (dT*dT / 2)*(p_frc[gg2]);
				


				//				FinRx << inV_vect[gg2].x / inV_vect[gg2].z << " ";
				//				FinRx << inR_vect[gg2].y << " ";
				//				FinRx << 0 << "\n";
				//
				//				FinVx << inV_vect[gg2].x << " ";
				//				FinVx << inV_vect[gg2].y << " ";
				//				FinVx << inV_vect[gg2].z << "\n";

				FinRx << outR_vect[gg2].x << " ";
				FinRx << outR_vect[gg2].y << " ";
				FinRx << outR_vect[gg2].z << "\n";

				//FinVx << outV_vect[gg2].x << " ";
		    	//FinVx << outV_vect[gg2].y << " ";
				//FinVx << outV_vect[gg2].z << "\n";
			}

		}
	}
	

	for (int gg22 = 0; gg22 < CP; gg22++)
	{
		//FoutRx << outR_vect[gg22].x << " ";
		//FoutRx << outR_vect[gg22].y << " ";
		//FoutRx << outR_vect[gg22].z << "\n";

		//FoutVx << outR_vect[gg22].x<< " ";
		//FoutVx << outV_vect[gg22].x / outV_vect[gg22].z << " ";		
		//FoutVx << 0 << "\n";

		FinVx << outR_vect[gg22].y << " ";
		
		FinVx << outV_vect[gg22].y / outV_vect[gg22].z << " ";	
		//FinVx << (outV_vect[gg22].y* outV_vect[gg22].z) / ((outV_vect[gg22].y * outV_vect[gg22].y) + (outV_vect[gg22].z * outV_vect[gg22].z)) << " ";
	    FinVx << 0 << "\n";

		//FoutVx << outV_vect[gg22].x << " ";
		//FoutVx << outV_vect[gg22].y << " ";		
		//FoutVx << outV_vect[gg22].z << "\n";
	}

	FinRx.close();
	FinVx.close();
	FoutRx.close();
	FoutVx.close();
	/*Gnuplot gp1;
	gp1.cmd( "set xrange[0:10]");
	gp1.cmd("set yrange[0:10]");
	gp1.cmd("set zrange[0:10]");
	gp1.cmd("splot 'inRx1.txt' with dots");

	//Gnuplot gp2;
	//gp2.cmd("set xrange[0:30]");
	//gp2.cmd("set yrange[0:30]");
	//gp2.cmd("set zrange[0:30]");
	//gp2.cmd("splot 'outRx1.txt' with dots");

	//Gnuplot gp3;
	//gp3.cmd("splot 'inVx1.txt' with dots");

	Gnuplot gp3;
	gp3.cmd("set xrange[-1:1]");
	gp3.cmd("set yrange[-1:1]");
	gp3.cmd("set zrange[-1:1]");
	gp3.cmd("splot 'inVx1.txt' with dots");

	Gnuplot gp4;
	gp4.cmd("set xrange[-1:1]");
	gp4.cmd("set yrange[-1:1]");
	gp4.cmd("set zrange[-1:1]");
	gp4.cmd("splot 'outVx1.txt' with dots");*/

	getchar();
	return 0;
}


