#include <iostream>
#include "SFML\Graphics.hpp"
#include "SFML\System\Time.hpp"
#include <vector>
#include <cassert>
#include <random>
#include "utilities.h"

#define PI 3.14159265

using namespace std;
using namespace sf;

#define E_M
//#define DIFF
//#define VEL_PROF
//#define VEL_DIST
#define GHOST_PART
#define GRID_SHIFT
#define POLYMER
//#define DRAW_PARTICLES
//#define SAVE_IMG 
//#define E2E_LEN

//Raw Values
const int M = 9;
const int a = 1;
const int Lx = 30 * a;		//Flow direction
const int Ly = 30 * a;
const int Lz = 30 * a;
const int m = 1;			//Solvent particle mass
const float kB_T = 0.011;
const float Tau = sqrt(m*a*a / kB_T);
const float delta_T = 0.15*Tau;
const float alpha = 130.0;
const float g = 0.000;

//Polymer properties
#ifdef POLYMER
const int n = 100;				//No. of monomers
const int m_n = M*m;			//Monomer mass
#else
const int n = 0;				//No. of monomers
const int m_n = 0;			//Monomer mass
#endif
const int deltaT_ratio = 100;
const float delta_T_n = (1.0/deltaT_ratio)*delta_T;
Color particleColor_n = Color(0, 191, 255);
const float sigma = a / 2.0;
const float sigma_2 = sigma*sigma;
const float epsilon = kB_T;
const float r_cut_2 = sigma*sigma*1.259921;
const float K = 500;		//Spring constant
const float l = a;		//Monomer equilibrium separation


//Derived values
const int N = M * Lx * Ly *Lz / (a*a*a);
const int N_tot = N + n;
const int dispX = 25;
const int dispY = 25;
const int WIDTH = dispX*Lx;
const int HEIGHT = dispX*Ly;
const float delta_T_2 = 0.5*delta_T*delta_T;
const float delta_T_2_n = 0.5*delta_T_n*delta_T_n;
const float radius = 1.0;
const float radius_n = sigma*dispX;
Color particleColor = Color(255, 255, 255);
const float rCos = cos(alpha*PI / 180);
const float rSin = sin(alpha*PI / 180);	
const float D = (kB_T*delta_T / (2 * m))*(-1 + (2 * M) / ((1-rCos)*(M-1+exp(-M))));

const float eta = kB_T*delta_T*(-0.5 + 5*M/((4-2*rCos-2*(2*rCos*rCos+1))*(M-1+exp(-M)))) +		//kinetic
					(m/(18*delta_T*M))*(M-1+exp(-M))*(1-rCos);									//collisional


mt19937 mt_rand(time(0));
normal_distribution<float> nd(0, 1);
uniform_real_distribution<float> ud(-1, 1);



float gd(float* _mean, float* _stddev)
{
	return ((*_mean) + (*_stddev) * nd(mt_rand));
}
float gd(float _mean, float _stddev)
{
	return ((_mean) + (_stddev) * nd(mt_rand));
}


struct Point2f
{
	float x;
	float y;

	Point2f operator + (Point2f a)
	{
		Point2f c;
		c.x = x + a.x;
		c.y = y + a.y;
		return c;
	}
};

class Particle
{
public:
	float sX;
	float sY;
	float sZ;
	float _sX;
	float _sY;
	float _sZ;

	int wX;
	int wY;
	int wZ;

	float vX;
	float vY;
	float vZ;

	float mass;

	int box_i;
	int box_j;
	int box_k;

	Particle()
	{
		mass = m;

		sX = Lx * (rand() / (float)RAND_MAX);
		sY = Ly * (rand() / (float)RAND_MAX);
		sZ = Lz * (rand() / (float)RAND_MAX);

		_sX = sX;
		_sY = sY;
		_sZ = sZ;

		wX = 0;
		wY = 0;
		wZ = 0;
		vX = (rand() / (float)RAND_MAX);
		vY = (rand() / (float)RAND_MAX);
		vZ = (rand() / (float)RAND_MAX);

		box_i = -1;
		box_j = -1;
		box_k = -1;
	}
};

class Box
{
public:
	float uX;
	float uY;
	float uZ;

	int count;			//Box solvent particles count
	int count_n;		//Box monomer count
	float Rx;
	float Ry;
	float Rz;

	double e_Kin;
	double scaleFactor;

	Box()
	{
		Reset();
	}

	void Reset()
	{
		uX = 0;
		uY = 0;
		uZ = 0;
		count = 0;
		count_n = 0;
		Rx = 0;
		Ry = 0;
		Rz = 0;
		e_Kin = 0;
		scaleFactor = 1;
	}
};

float LJPa(float _rs)
{
	if (_rs > r_cut_2)
		return 0.0;
	else
	{
		//Assuming epsilon or a = 1
		//float _r = sqrt(_rs);
		return (epsilon / m_n)*((48 * pow(sigma_2 / _rs, 6) / _rs) - (24 * pow(sigma_2 / _rs, 3) / _rs));
	}
}

void Initialize(Particle particles[])
{
#ifdef POLYMER
	//Position initialization for polymer
	particles[0].mass = m_n;
	particles[0].sX = 2.5;
	particles[0].sY = Ly / 2;
	particles[0].sZ = Lz / 2;
	particles[0]._sX = particles[0].sX;
	particles[0]._sY = particles[0].sY;
	particles[0]._sZ = particles[0].sZ;

	int stepX = 1 * (2 * sigma);
	int prev_stepX = stepX;
	int stepY = 0;
	for (int i = 1; i < n; i++)
	{
		particles[i].mass = m_n;
		if ((particles[i - 1].sX >= Lx - 2.5 && stepX > 0) || (particles[i - 1].sX <= 2.5 && stepX < 0))
		{
			prev_stepX = stepX;
			stepX = 0;
			stepY += 1;

		}

		particles[i].sX = particles[i - 1].sX + stepX;
		particles[i].sY = Ly / 2 + stepY;
		particles[i].sZ = Lz / 2;

		particles[i]._sX = particles[i].sX;
		particles[i]._sY = particles[i].sY;
		particles[i]._sZ = particles[i].sZ;

		if (stepX == 0)
		{
			stepX = -prev_stepX;
			prev_stepX = stepX;
		}
	}
#endif

	//Other initializations
	float uX = 0;
	float uY = 0;
	float uZ = 0;

	for (int i = 0; i < N_tot; i++)
	{
		uX += particles[i].vX*particles[i].mass;
		uY += particles[i].vY*particles[i].mass;
		uZ += particles[i].vZ*particles[i].mass;
	}

	//Assuming equal masses
	uX /= N*m + n*m_n;
	uY /= N*m + n*m_n;
	uZ /= N*m + n*m_n;

	//Velocity Normalization
	for (int i = 0; i < N_tot; i++)
	{
		particles[i].vX -= uX;
		particles[i].vY -= uY;
		particles[i].vZ -= uZ;
	}

	//Temperature Normalization
	float curr_kB_T = 0;
	for (int i = 0; i < N_tot; i++)
	{
		curr_kB_T += particles[i].mass*(particles[i].vX*particles[i].vX + particles[i].vY*particles[i].vY + 
							particles[i].vZ*particles[i].vZ);
	}
	curr_kB_T /= 2 * N_tot;
	
	if (curr_kB_T == 0)
	{
		cout << "Error currTemp = 0 !" << endl;
		return;
	}

	float ratio = sqrt(kB_T / curr_kB_T);
	for (int i = 0; i < N_tot; i++)
	{
		particles[i].vX *= ratio;
		particles[i].vY *= ratio;
		particles[i].vZ *= ratio;
	}
	
	//Add desired velocity to system ...
	/*float add_vX = 2;
	float add_vY = 0;
	for (int i = 0; i < N; i++)
	{
		particles[i].vX += add_vX;
		particles[i].vY += add_vY;
	}*/
	//Setting CircleShapes class
}

void BC_Periodic(Particle& particle)
{
	if (particle.sX >= Lx || particle.sX < 0)
	{
		int warpX = floor(particle.sX / Lx);
		particle.wX += warpX;
		particle.sX -= warpX * Lx;
	}
	if (particle.sY >= Ly || particle.sY < 0)
	{
		int warpY = floor(particle.sY / Ly);
		particle.wY += warpY;
		particle.sY -= warpY * Ly;
	}
	if (particle.sZ >= Lz || particle.sZ < 0)
	{
		int warpZ = floor(particle.sZ / Lz);
		particle.wZ += warpZ;
		particle.sZ -= warpZ * Lz;
	}
}

void BC_Slip(Particle& particle)
{
	//Assume after one reflection particle stays within bounds
	if (particle.sX >= Lx)
	{
		particle.sX -= 2 * (particle.sX - Lx);
		particle.vX *= -1;
	}
	else if (particle.sX < 0)
	{
		particle.sX -= 2 * particle.sX;
		particle.vX *= -1;
	}
	
	if (particle.sY >= Ly)
	{
		particle.sY -= 2 * (particle.sY - Ly);
		particle.vY *= -1;
	}
	else if (particle.sY < 0)
	{
		particle.sY -= 2 * particle.sY;
		particle.vY *= -1;
	}

	if (particle.sZ >= Lz)
	{
		particle.sZ -= 2 * (particle.sZ - Lz);
		particle.vZ *= -1;
	}
	else if (particle.sZ < 0)
	{
		particle.sZ -= 2 * particle.sZ;
		particle.vZ *= -1;
	}
}

void BC_Poiseuille(Particle& particle)
{
	// Z - no slip
	if (particle.sZ >= Lz)
	{
		particle.sZ -= 2 * (particle.sZ - Lz);
		float _time = (Lz - particle._sZ) / particle.vZ;
		
		particle.sX = particle._sX + particle.vX * _time + 0.5*g*_time*_time;
		particle.vX *= -1;
		particle.sX = particle.sX + particle.vX * (delta_T - _time) + 0.5*g*(delta_T - _time)*(delta_T - _time);
		particle.sY = particle._sY + particle.vY*(2 * _time - delta_T);
		particle.vY *= -1;
		particle.vZ *= -1;

	}
	else if (particle.sZ <= 0)
	{
		particle.sZ -= 2 * particle.sZ;
		float _time = (-particle._sZ) / particle.vZ;
		
		particle.sX = particle._sX + particle.vX * _time + 0.5*g*_time*_time;
		particle.vX *= -1;
		particle.sX = particle.sX + particle.vX * (delta_T - _time) + 0.5*g*(delta_T - _time)*(delta_T - _time);
		particle.sY = particle._sY + particle.vY*(2 * _time - delta_T);
		particle.vY *= -1;
		particle.vZ *= -1;
	}

	// X - periodic
	if (particle.sX >= Lx || particle.sX < 0)
	{
		int warpX = floor(particle.sX / Lx);
		particle.wX += warpX;
		particle.sX -= warpX * Lx;
	}

	// Y - periodic
	if (particle.sY >= Ly || particle.sY < 0)
	{
		int warpY = floor(particle.sY / Ly);
		particle.wY += warpY;
		particle.sY -= warpY * Ly;
	}
}

void Stream(Particle particles[])
{
	//Solvent ballistic streaming step
	for (int i = n; i < N_tot; i++)
	{
		particles[i]._sX = particles[i].sX;
		particles[i]._sY = particles[i].sY;
		particles[i]._sZ = particles[i].sZ;
		particles[i].sX += particles[i].vX * delta_T + g*delta_T_2;
		particles[i].sY += particles[i].vY * delta_T;
		particles[i].sZ += particles[i].vZ * delta_T;

		//BC_Periodic(particles[i]);
		//BC_Slip(particles[i]);
		BC_Poiseuille(particles[i]);

		particles[i].vX += g*delta_T;
	}

#ifdef POLYMER
	//Molecular dynamics step
	float AccX[n][n] = { 0 };
	float AccY[n][n] = { 0 };
	float AccZ[n][n] = { 0 };
	float AccX_[n][n] = { 0 };
	float AccY_[n][n] = { 0 };
	float AccZ_[n][n] = { 0 };
	float r_K[n - 1] = { 0 };
	float K_Acc_X[n+1] = { 0 };
	float K_Acc_Y[n+1] = { 0 };
	float K_Acc_Z[n+1] = { 0 };
	float K_Acc_X_[n+1] = { 0 };
	float K_Acc_Y_[n+1] = { 0 };
	float K_Acc_Z_[n+1] = { 0 };
	float x_ = 0;
	float y_ = 0;
	float z_ = 0;

	//Initial Force due to Lennard-Jonnes
	for (int i = 0; i < n; i++)
	{
		for (int j = i + 2; j < n; j++)
		{
			if (abs(particles[i].sX - particles[j].sX) <= Lx / 2.0)			//sX
				x_ = particles[i].sX - particles[j].sX;
			else
			{
				if (particles[i].sX >= particles[j].sX)
					x_ = -Lx + (particles[i].sX - particles[j].sX);
				else
					x_ = Lx + (particles[i].sX - particles[j].sX);
			}

			if (abs(particles[i].sY - particles[j].sY) <= Ly / 2.0)			//sY
				y_ = particles[i].sY - particles[j].sY;
			else
			{
				if (particles[i].sY >= particles[j].sY)
					y_ = -Ly + (particles[i].sY - particles[j].sY);
				else
					y_ = Ly + (particles[i].sY - particles[j].sY);
			}

			if (abs(particles[i].sZ - particles[j].sZ) <= Lz / 2.0)			//sZ
				z_ = particles[i].sZ - particles[j].sZ;
			else
			{
				if (particles[i].sZ >= particles[j].sZ)
					z_ = -Lz + (particles[i].sZ - particles[j].sZ);
				else
					z_ = Lz + (particles[i].sZ - particles[j].sZ);
			}

			float dist = x_*x_ + y_*y_ + z_*z_;
			float acc = LJPa(dist);
			AccX[i][j] = acc*x_;
			AccX[j][i] = -AccX[i][j];
			AccY[i][j] = acc*y_;
			AccY[j][i] = -AccY[i][j];
			AccZ[i][j] = acc*z_;
			AccZ[j][i] = -AccZ[i][j];
		}
	}
	//Initial Force due to spring
	for (int i = 0; i < n-1; i++)
	{
		if (abs(particles[i].sX - particles[i + 1].sX) <= Lx / 2.0)			//sX
			x_ = particles[i + 1].sX - particles[i].sX;
		else
		{
			if (particles[i].sX >= particles[i + 1].sX)
				x_ = Lx + (particles[i + 1].sX - particles[i].sX);
			else
				x_ = -Lx + (particles[i + 1].sX - particles[i].sX);
		}

		if (abs(particles[i].sY - particles[i + 1].sY) <= Ly / 2.0)			//sY
			y_ = particles[i + 1].sY - particles[i].sY;
		else
		{
			if (particles[i].sY >= particles[i + 1].sY)
				y_ = Ly + (particles[i + 1].sY - particles[i].sY);
			else
				y_ = -Ly + (particles[i + 1].sY - particles[i].sY);
		}

		if (abs(particles[i].sZ - particles[i + 1].sZ) <= Lz / 2.0)			//sZ
			z_ = particles[i + 1].sZ - particles[i].sZ;
		else
		{
			if (particles[i].sZ >= particles[i + 1].sZ)
				z_ = Lz + (particles[i + 1].sZ - particles[i].sZ);
			else
				z_ = -Lz + (particles[i + 1].sZ - particles[i].sZ);
		}
		
		r_K[i] = sqrt(x_*x_ + y_*y_ + z_*z_);
		K_Acc_X[i + 1] = (K / m_n)*(1 - l / r_K[i])*x_;
		K_Acc_Y[i + 1] = (K / m_n)*(1 - l / r_K[i])*y_;
		K_Acc_Z[i + 1] = (K / m_n)*(1 - l / r_K[i])*z_;
	}
	
	// MD loop
	for (int k = 0; k < deltaT_ratio; k++)
	{
		//Update positions
		for (int i = 0; i < n; i++)
		{
			particles[i]._sX = particles[i].sX;
			particles[i]._sY = particles[i].sY;
			particles[i]._sZ = particles[i].sZ;
			float accX = 0;
			float accY = 0;
			float accZ = 0;
			for (int j = 0; j < n; j++)
			{
				accX += AccX[i][j];
				accY += AccY[i][j];
				accZ += AccZ[i][j];
			}
			particles[i].sX += particles[i].vX * delta_T_n + (g + accX + K_Acc_X[i + 1] - K_Acc_X[i])*delta_T_2_n;
			particles[i].sY += particles[i].vY * delta_T_n + (accY + K_Acc_Y[i + 1] - K_Acc_Y[i])*delta_T_2_n;
			particles[i].sZ += particles[i].vZ * delta_T_n + (accZ + K_Acc_Z[i + 1] - K_Acc_Z[i])*delta_T_2_n;
			
			BC_Periodic(particles[i]);
		}
		
		//Find new accelerations
		// Force due to Lennard-Jones
		for (int i = 0; i < n; i++)
		{
			for (int j = i + 2; j < n; j++)
			{
				if (abs(particles[i].sX - particles[j].sX) <= Lx / 2.0)			//sX
					x_ = particles[i].sX - particles[j].sX;
				else
				{
					if (particles[i].sX >= particles[j].sX)
						x_ = -Lx + (particles[i].sX - particles[j].sX);
					else
						x_ = Lx + (particles[i].sX - particles[j].sX);
				}

				if (abs(particles[i].sY - particles[j].sY) <= Ly / 2.0)			//sY
					y_ = particles[i].sY - particles[j].sY;
				else
				{
					if (particles[i].sY >= particles[j].sY)
						y_ = -Ly + (particles[i].sY - particles[j].sY);
					else
						y_ = Ly + (particles[i].sY - particles[j].sY);
				}

				if (abs(particles[i].sZ - particles[j].sZ) <= Lz / 2.0)			//sZ
					z_ = particles[i].sZ - particles[j].sZ;
				else
				{
					if (particles[i].sZ >= particles[j].sZ)
						z_ = -Lz + (particles[i].sZ - particles[j].sZ);
					else
						z_ = Lz + (particles[i].sZ - particles[j].sZ);
				}

				float dist = x_*x_ + y_*y_ + z_*z_;
				float acc = LJPa(dist);
				AccX_[i][j] = acc*x_;
				AccX_[j][i] = -AccX_[i][j];
				AccY_[i][j] = acc*y_;
				AccY_[j][i] = -AccY_[i][j];
				AccZ_[i][j] = acc*z_;
				AccZ_[j][i] = -AccZ_[i][j];
			}
		}

		// Force due to Spring
		for (int i = 0; i < n - 1; i++)
		{
			if (abs(particles[i].sX - particles[i + 1].sX) <= Lx / 2.0)			//sX
				x_ = particles[i + 1].sX - particles[i].sX;
			else
			{
				if (particles[i].sX >= particles[i + 1].sX)
					x_ = Lx + (particles[i + 1].sX - particles[i].sX);
				else
					x_ = -Lx + (particles[i + 1].sX - particles[i].sX);
			}

			if (abs(particles[i].sY - particles[i + 1].sY) <= Ly / 2.0)			//sY
				y_ = particles[i + 1].sY - particles[i].sY;
			else
			{
				if (particles[i].sY >= particles[i + 1].sY)
					y_ = Ly + (particles[i + 1].sY - particles[i].sY);
				else
					y_ = -Ly + (particles[i + 1].sY - particles[i].sY);
			}

			if (abs(particles[i].sZ - particles[i + 1].sZ) <= Lz / 2.0)			//sZ
				z_ = particles[i + 1].sZ - particles[i].sZ;
			else
			{
				if (particles[i].sZ >= particles[i + 1].sZ)
					z_ = Lz + (particles[i + 1].sZ - particles[i].sZ);
				else
					z_ = -Lz + (particles[i + 1].sZ - particles[i].sZ);
			}

			r_K[i] = sqrt(x_*x_ + y_*y_ + z_*z_);
			K_Acc_X_[i + 1] = (K / m_n)*(1 - l / r_K[i])*x_;
			K_Acc_Y_[i + 1] = (K / m_n)*(1 - l / r_K[i])*y_;
			K_Acc_Z_[i + 1] = (K / m_n)*(1 - l / r_K[i])*z_;
		}
		

		//Update velocities
		for (int i = 0; i < n; i++)
		{
			float accX = 0;
			float accY = 0;
			float accZ = 0;
			float accX_ = 0;
			float accY_ = 0;
			float accZ_ = 0;
			for (int j = 0; j < n; j++)
			{
				accX += AccX[i][j];
				accX_ += AccX_[i][j];
				accY += AccY[i][j];
				accY_ += AccY_[i][j];
				accZ += AccZ[i][j];
				accZ_ += AccZ_[i][j];
			}

			particles[i].vX += (g + 0.5*(accX + accX_ + (K_Acc_X[i + 1] - K_Acc_X[i]) + (K_Acc_X_[i + 1] - K_Acc_X_[i])))*delta_T_n;
			particles[i].vY += (0.5*(accY + accY_ + (K_Acc_Y[i + 1] - K_Acc_Y[i]) + (K_Acc_Y_[i + 1] - K_Acc_Y_[i])))*delta_T_n;
			particles[i].vZ += (0.5*(accZ + accZ_ + (K_Acc_Z[i + 1] - K_Acc_Z[i]) + (K_Acc_Z_[i + 1] - K_Acc_Z_[i])))*delta_T_n;
		}
		//Copy accelerations
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				AccX[i][j] = AccX_[i][j];
				AccY[i][j] = AccY_[i][j];
				AccZ[i][j] = AccZ_[i][j];
			}

			K_Acc_X[i+1] = K_Acc_X_[i+1];
			K_Acc_Y[i+1] = K_Acc_Y_[i+1];
			K_Acc_Z[i+1] = K_Acc_Z_[i+1];
		}
	}
#endif

	return;
}

void SRDcollide(Particle particles[], Box boxes[][Ly + 2][Lx + 2])
{
	for (int i = 0; i < Lz + 2; i++)
	{
		for (int j = 0; j < Ly + 2; j++)
		{
			for (int k = 0; k < Lx + 2; k++)
			{
				boxes[i][j][k].Reset();
			}
		}
	}

#ifdef GRID_SHIFT
	//Shift the grid
	float shiftX = ud(mt_rand)*a / 2;
	float shiftY = ud(mt_rand)*a / 2;
	float shiftZ = ud(mt_rand)*a / 2;
#else
	float shiftX = 0;
	float shiftY = 0;
	float shiftZ = 0;
#endif

	int box_Z = 0;	// z-box cutting the wall
	if (shiftZ > 0)
		box_Z = Lz;
	else
		box_Z = 1;


	for (int i = 0; i < N_tot; i++)
	{
		int K = (int)(particles[i].sX + shiftX + 1);
		int J = (int)(particles[i].sY + shiftY + 1);
		int I = (int)(particles[i].sZ + shiftZ + 1);
		particles[i].box_i = I;
		particles[i].box_j = J;
		particles[i].box_k = K;

		if (i < n)
			boxes[I][J][K].count_n++;
		else
			boxes[I][J][K].count++;

		boxes[I][J][K].uX += particles[i].vX*particles[i].mass;
		boxes[I][J][K].uY += particles[i].vY*particles[i].mass;
		boxes[I][J][K].uZ += particles[i].vZ*particles[i].mass;
	}
	for (int i = 0; i < Lz + 2; i++)
	{
		for (int j = 0; j < Ly + 2; j++)
		{
			for (int k = 0; k < Lx + 2; k++)
			{
				if (boxes[i][j][k].count + boxes[i][j][k].count_n != 0)
				{
#ifdef GHOST_PART

					if ((boxes[i][j][k].count < M) && (i == 0 || i == Lz + 1 || i == box_Z))
					{
						float stddev = sqrt((M - boxes[i][j][k].count)*kB_T / m);
						boxes[i][j][k].uX = (boxes[i][j][k].uX + m*gd(0, stddev)) / (M*m + boxes[i][j][k].count_n*m_n);
						boxes[i][j][k].uY = (boxes[i][j][k].uY + m*gd(0, stddev)) / (M*m + boxes[i][j][k].count_n*m_n);
						boxes[i][j][k].uZ = (boxes[i][j][k].uZ + m*gd(0, stddev)) / (M*m + boxes[i][j][k].count_n*m_n);
						boxes[i][j][k].Rx = ud(mt_rand);
						boxes[i][j][k].Ry = ud(mt_rand);
						boxes[i][j][k].Rz = ud(mt_rand);
						float mag = sqrt(boxes[i][j][k].Rx*boxes[i][j][k].Rx +
										 boxes[i][j][k].Ry*boxes[i][j][k].Ry +
										 boxes[i][j][k].Rz*boxes[i][j][k].Rz);
						boxes[i][j][k].Rx /= mag;
						boxes[i][j][k].Ry /= mag;
						boxes[i][j][k].Rz /= mag;
					}
					else
					{
						boxes[i][j][k].uX /= m*boxes[i][j][k].count + m_n*boxes[i][j][k].count_n;
						boxes[i][j][k].uY /= m*boxes[i][j][k].count + m_n*boxes[i][j][k].count_n;
						boxes[i][j][k].uZ /= m*boxes[i][j][k].count + m_n*boxes[i][j][k].count_n;
						boxes[i][j][k].Rx = ud(mt_rand);
						boxes[i][j][k].Ry = ud(mt_rand);
						boxes[i][j][k].Rz = ud(mt_rand);
						float mag = sqrt(boxes[i][j][k].Rx*boxes[i][j][k].Rx +
							boxes[i][j][k].Ry*boxes[i][j][k].Ry +
							boxes[i][j][k].Rz*boxes[i][j][k].Rz);
						boxes[i][j][k].Rx /= mag;
						boxes[i][j][k].Ry /= mag;
						boxes[i][j][k].Rz /= mag;
					}
#else
					boxes[i][j][k].uX /= m*boxes[i][j][k].count + m_n*boxes[i][j][k].count_n;
					boxes[i][j][k].uY /= m*boxes[i][j][k].count + m_n*boxes[i][j][k].count_n;
					boxes[i][j][k].uZ /= m*boxes[i][j][k].count + m_n*boxes[i][j][k].count_n;
					boxes[i][j][k].Rx = ud(mt_rand);
					boxes[i][j][k].Ry = ud(mt_rand);
					boxes[i][j][k].Rz = ud(mt_rand);
					float mag = sqrt(boxes[i][j][k].Rx*boxes[i][j][k].Rx +
						boxes[i][j][k].Ry*boxes[i][j][k].Ry +
						boxes[i][j][k].Rz*boxes[i][j][k].Rz);
					boxes[i][j][k].Rx /= mag;
					boxes[i][j][k].Ry /= mag;
					boxes[i][j][k].Rz /= mag;
#endif
				}
			}
		}
	}


	//Thermalisation incorporated in the collision step...
	/*for (int i = 0; i < N; i++)
	{
		int I = particles[i].box_i;
		int J = particles[i].box_j;

		double _dVx = particles[i].vX - boxes[I][J].uX;
		double _dVy = particles[i].vY - boxes[I][J].uY;
		boxes[I][J].e_Kin += 0.5*m*(_dVx*_dVx + _dVy*_dVy);
	}
	for (int i = 1; i < Ly + 1; i++)
	{
		for (int j = 1; j < Lx + 1; j++)
		{
			if (boxes[i][j].count > 1 && boxes[i][j].e_Kin != 0)
			{
				boxes[i][j].scaleFactor = sqrt(kB_T * 1 * (boxes[i][j].count - 1) / boxes[i][j].e_Kin);

				boxes[i][j].scaleFactor = sqrt(gd(kB_T*(boxes[i][j].count - 1), kB_T*sqrt(boxes[i][j].count - 1)) / boxes[i][j].e_Kin);
				if (abs(boxes[i][j].scaleFactor) > 100)
					cout << boxes[i][j].e_Kin << endl;
			}
		}
	}*/

	//Actual collision step here ...
	for (int i = 0; i < N_tot; i++)
	{
		int I = particles[i].box_i;
		int J = particles[i].box_j;
		int K = particles[i].box_k;

		float dvX = particles[i].vX - boxes[I][J][K].uX;
		float dvY = particles[i].vY - boxes[I][J][K].uY;
		float dvZ = particles[i].vZ - boxes[I][J][K].uZ;
		
		float dot = dvX*boxes[I][J][K].Rx + dvY*boxes[I][J][K].Ry + dvZ*boxes[I][J][K].Rz;
		float dvX_ll = boxes[I][J][K].Rx*dot;
		float dvY_ll = boxes[I][J][K].Ry*dot;
		float dvZ_ll = boxes[I][J][K].Rz*dot;
		float dvX_p = dvX - dvX_ll;
		float dvY_p = dvY - dvY_ll;
		float dvZ_p = dvZ - dvZ_ll;

		particles[i].vX = boxes[I][J][K].uX + dvX_ll + dvX_p*rCos + 
						 (dvY_p*boxes[I][J][K].Rz - dvZ_p*boxes[I][J][K].Ry)*rSin;
		particles[i].vY = boxes[I][J][K].uY + dvY_ll + dvY_p*rCos +
						 (dvZ_p*boxes[I][J][K].Rx - dvX_p*boxes[I][J][K].Rz)*rSin;
		particles[i].vZ = boxes[I][J][K].uZ + dvZ_ll + dvZ_p*rCos +
						 (dvX_p*boxes[I][J][K].Ry - dvY_p*boxes[I][J][K].Rx)*rSin;
	}

	return;
}

float Energy(Particle particles[], int begin = 0, int end = N_tot)
{
	double energy = 0;
	for (int i = begin; i < end; i++)
	{
		energy += 0.5*particles[i].mass*
				  (particles[i].vX*particles[i].vX + 
				   particles[i].vY*particles[i].vY + 
				   particles[i].vZ*particles[i].vZ);
	}
	
	return energy/(end-begin);
}

float Momentum(Particle particles[], int begin = 0, int end = N_tot)
{
	double momentumX = 0;
	double momentumY = 0;
	double momentumZ = 0;
	for (int i = begin; i < end; i++)
	{
		momentumX += particles[i].mass*particles[i].vX;
		momentumY += particles[i].mass*particles[i].vY;
		momentumZ += particles[i].mass*particles[i].vZ;
	}
	/*momentumX = momentumX / N;
	momentumY = momentumY / N;
	momentumZ = momentumZ / N;*/
	//cout << sqrt(momentumX*momentumX + momentumY*momentumY) << endl;

	return sqrt(momentumX*momentumX + momentumY*momentumY + momentumZ*momentumZ);
}

void Thermostat(Particle particles[])		// !!!!!
{
	double e_Kin = 0.0;
	double vX_Cm = 0.0;
	double vY_Cm = 0.0;
	double vZ_Cm = 0.0;
	for (int i = 0; i < N; i++)
	{
		vX_Cm += particles[i].vX;
		vY_Cm += particles[i].vY;
		vZ_Cm += particles[i].vZ;
	}
	vX_Cm /= N;
	vY_Cm /= N;
	vZ_Cm /= N;
	for (int i = 0; i < N; i++)
	{
		e_Kin += 0.5 * ((particles[i].vX - vX_Cm)*(particles[i].vX - vX_Cm) +
						(particles[i].vY - vY_Cm)*(particles[i].vY - vY_Cm) +
						(particles[i].vZ - vZ_Cm)*(particles[i].vZ - vZ_Cm));
	}
	double therFactor = sqrt(kB_T * 3 * (N - 1) / (2 * e_Kin));
	for (int i = 0; i < N; i++)
	{
		particles[i].vX *= therFactor;
		particles[i].vY *= therFactor;
		particles[i].vZ *= therFactor;
	}
	return;
}

void Thermostat_Local(Particle particles[], Box boxes[][Ly + 2][Lx + 2])
{
	for (int i = 0; i < Lz + 2; i++)
	{
		for (int j = 0; j < Ly + 2; j++)
		{
			for (int k = 0; k < Lx + 2; k++)
			{
				boxes[i][j][k].Reset();
			}
		}
	}

	for (int i = n; i < N_tot; i++)
	{
		int K = (int)(particles[i].sX + 1);
		int J = (int)(particles[i].sY + 1);
		int I = (int)(particles[i].sZ + 1);
		particles[i].box_i = I;
		particles[i].box_j = J;
		particles[i].box_k = K;

		/*if (i < n)
			boxes[I][J][K].count_n++;
		else
			boxes[I][J][K].count++;*/

		boxes[I][J][K].uX += particles[i].vX*particles[i].mass;
		boxes[I][J][K].uY += particles[i].vY*particles[i].mass;
		boxes[I][J][K].uZ += particles[i].vZ*particles[i].mass;
	}

	for (int i = 0; i < Lz + 2; i++)
	{
		for (int j = 0; j < Ly + 2; j++)
		{
			for (int k = 0; k < Lx + 2; k++)
			{
				if (boxes[i][j][k].count + boxes[i][j][k].count_n != 0)
				{
					float m_effective = boxes[i][j][k].count*m + boxes[i][j][k].count_n*m_n;
					boxes[i][j][k].uX /= m_effective;
					boxes[i][j][k].uY /= m_effective;
					boxes[i][j][k].uZ /= m_effective;
				}

			}
		}
	}

	//Assuming equal masses
	for (int i = n; i < N_tot; i++)
	{
		double _dVx = particles[i].vX - boxes[particles[i].box_i][particles[i].box_j][particles[i].box_k].uX;
		double _dVy = particles[i].vY - boxes[particles[i].box_i][particles[i].box_j][particles[i].box_k].uY;
		double _dVz = particles[i].vZ - boxes[particles[i].box_i][particles[i].box_j][particles[i].box_k].uZ;
		boxes[particles[i].box_i][particles[i].box_j][particles[i].box_k].e_Kin += 
												0.5*particles[i].mass*(_dVx*_dVx + _dVy*_dVy + _dVz*_dVz);
	}

	for (int i = 1; i < Lz + 1; i++)
	{
		for (int j = 1; j < Ly + 1; j++)
		{
			for (int k = 1; k < Lx + 1; k++)
			{
				if (boxes[i][j][k].count + boxes[i][j][k].count_n > 1 && boxes[i][j][k].e_Kin != 0)
				{
					boxes[i][j][k].scaleFactor = sqrt(kB_T * 3 * (boxes[i][j][k].count + boxes[i][j][k].count_n - 1) /
						(2 * boxes[i][j][k].e_Kin));

					boxes[i][j][k].scaleFactor =
						sqrt(abs(gd(kB_T * 3 * (boxes[i][j][k].count + boxes[i][j][k].count_n - 1), 
						kB_T*sqrt(boxes[i][j][k].count + boxes[i][j][k].count_n - 1))) /
						(2 * boxes[i][j][k].e_Kin));
				}
			}
		}
	}

	for (int i = n; i < N_tot; i++)
	{
		int I = particles[i].box_i;
		int J = particles[i].box_j;
		int K = particles[i].box_k;

		if (boxes[I][J][K].count + boxes[I][J][K].count_n > 1 && boxes[I][J][K].e_Kin != 0)
		{	
			particles[i].vX = boxes[I][J][K].uX + (particles[i].vX - boxes[I][J][K].uX)*boxes[I][J][K].scaleFactor;
			particles[i].vY = boxes[I][J][K].uY + (particles[i].vY - boxes[I][J][K].uY)*boxes[I][J][K].scaleFactor;
			particles[i].vZ = boxes[I][J][K].uZ + (particles[i].vZ - boxes[I][J][K].uZ)*boxes[I][J][K].scaleFactor;
		}
	}
	return;
}

void d_ParticlesCount(Box boxes[][Ly + 2][Lx + 2])
{
	int _sum = 0;
	for (int i = 1; i < Lz + 1; i++)
	{
		for (int j = 1; j < Ly + 1; j++)
		{
			for (int k = 1; k < Lx + 1; k++)
			{
				_sum += boxes[i][j][k].count;
			}
		}
	}
	
	cout << "Total number of particles : " << _sum << endl;
	return;
}


int main(int argc, char** argv[])
{
	assert(Lx % a == 0);		//Lx is integral multiple of a
	assert(Ly % a == 0);		//Ly is integral multiple of a
	assert(Lz % a == 0);		
	srand(time(NULL));
	

	cout << "Viscosity theoritical: " << eta << "\n\n";
	//cout << "Time length unit (Tau) : " << Tau << "\n\n";
	//cout << "Diffusion Coefficient : " << D << "\n\n";
	//system("pause");


	Particle particles[N_tot];
	Initialize(particles);

	Box boxes[Lz + 2][Ly + 2][Lx + 2];

	ContextSettings settings;
	settings.antialiasingLevel = 0;
	RenderWindow renderWindow(VideoMode(WIDTH, HEIGHT), "MPC");
	renderWindow.setPosition(Vector2i(20, 50));
	RenderWindow plotWindow(VideoMode(WIDTH, HEIGHT), "Statistics");
	plotWindow.setPosition(Vector2i(700, 50));
	renderWindow.setFramerateLimit(100);
	renderWindow.requestFocus();
	Clock fpsClock;
	
	RectangleShape lineX(Vector2f(HEIGHT, 1));
	lineX.setPosition(20, 0);
	lineX.rotate(90);
	RectangleShape lineY(Vector2f(WIDTH, 1));
	lineY.setPosition(0, HEIGHT - 100);
	lineY.rotate(0);
	const int hist_size = 100;
	int display_mode = 3;		
	int iter = -1;

#ifdef DRAW_PARTICLES
	CircleShape SolventParticles;
	SolventParticles.setRadius(radius);
	SolventParticles.setFillColor(particleColor);

	CircleShape MonomerParticles;
	MonomerParticles.setRadius(radius_n*0.9);
	MonomerParticles.setFillColor(particleColor_n);
	MonomerParticles.setOutlineThickness(radius_n*0.1);
	MonomerParticles.setOutlineColor(Color(255, 255, 255));

#ifdef SAVE_IMG
	const int anim_size = 2000;
	const int anim_start = 0;
	Vector2u windowSize = renderWindow.getSize();
	Texture texture[anim_size];
	for (int i = 0; i < anim_size; i++)
		texture[i].create(windowSize.x, windowSize.y);

	char outfile[20] = { ' ' };
#endif

#endif

#ifdef DIFF
	//For Diffusion Measurement
	const int diff_iters = 1000;
	const int diff_part = N_tot;
	assert(diff_part <= N_tot);
	float init_X[diff_part] = { 0 };
	float diff_X[diff_iters] = { 0 };
	float init_Y[diff_part] = { 0 };
	float diff_Y[diff_iters] = { 0 };
	float init_Z[diff_part] = { 0 };
	float diff_Z[diff_iters] = { 0 };
	for (int i = 0; i < diff_part; i++)
	{
		init_X[i] = particles[i].sX;
		init_Y[i] = particles[i].sY;
		init_Z[i] = particles[i].sZ;
	}
#endif

#ifdef VEL_PROF
	//For Velocity Profile - Poiseuille
	long double velProfileX_Z[Lz] = { 0 };
	int velProfileCountX_Z[Lz] = { 0 };
	long double velProfileX_Y[Ly] = { 0 };
	int velProfileCountX_Y[Ly] = { 0 };
	CircleShape velPoint(1);
	velPoint.setFillColor(Color(100, 255, 100));
	const int _maxIters = 1500;
	const int _y = 10;	//Y-Layer number
	assert(_y < Ly);
#endif

#ifdef E_M
	//For Energy Momentum Plots
	const int EM_max_iter = 500;
	const int EM_offset_iter = 50;
	float energyArr[EM_max_iter] = { 0 };
	float momentumArr[EM_max_iter] = { 0 };
#endif

#if defined(E2E_LEN) && defined(POLYMER)

	const int e2e_max_iter = 800;
	float e2e_len[e2e_max_iter] = { 0 };

#endif

///////////////////////////////////////////////////////////////////
	while (renderWindow.isOpen())
	{
		iter++;
		//int fps = 1.0 / (fpsClock.restart()).asSeconds();
		//cout << fps << endl;


		Stream(particles);
		SRDcollide(particles, boxes);
		//Thermostat(particles);
		Thermostat_Local(particles, boxes);


#ifdef DRAW_PARTICLES
		Event event;
		while (renderWindow.pollEvent(event))
		{
			if (event.type == Event::Closed)
				renderWindow.close();
			else if (event.type == Event::KeyPressed)
			{
				switch (event.key.code)
				{
				case Keyboard::Num1:display_mode = 1; renderWindow.create(VideoMode(dispX*Ly, dispX*Lz), "MPC"); break;
				case Keyboard::Num2:display_mode = 2; renderWindow.create(VideoMode(dispX*Lz, dispX*Lx), "MPC"); break;
				case Keyboard::Num3:display_mode = 3; renderWindow.create(VideoMode(dispX*Lx, dispX*Ly), "MPC"); break;
				default: break;
				}
			}
		}

		//Draw the particles
		renderWindow.clear();
		for (int i = 0; i < n; i++)
		{
			switch (display_mode)
			{
			case 1:MonomerParticles.setPosition(particles[i].sY * dispX, particles[i].sZ * dispX); break;
			case 2:MonomerParticles.setPosition(particles[i].sZ * dispX, particles[i].sX * dispX); break;
			case 3:MonomerParticles.setPosition(particles[i].sX * dispX, particles[i].sY * dispX); break;
			default:MonomerParticles.setPosition(particles[i].sY * dispX, particles[i].sZ * dispX); break;
			}
			renderWindow.draw(MonomerParticles);

		}

#ifdef SAVE_IMG
		if (iter < anim_size + anim_start && iter >= anim_start)
			texture[iter - 0].update(renderWindow);
		if (iter == anim_size + anim_start)
		{
			cout << "Saving captured images ... \n\n";
			for (int img = 0; img < anim_size; img++)
			{
				sprintf(outfile, "anim/shot_%d.png", img+1);
				texture[img].copyToImage().saveToFile(outfile);
			}
			cout << "Saving complete !!! \n\n";
		}
#endif

		/*for (int i = n; i < N_tot; i++)
		{
			SolventParticles.setPosition(particles[i].sX * dispX, particles[i].sY * dispY);
			renderWindow.draw(SolventParticles);
		}*/
		renderWindow.display();
#endif

		
#ifdef VEL_DIST
		/*Vertex line[] =
		{
			Vertex(Vector2f(WIDTH, 10)),
			Vertex(Vector2f(150, 150))
		};
		plotWindow.draw(line, 2, Lines);*/
		if (iter % 10 == 0)
		{
			float speed[N] = { 0 };
			float max_F = 0;
			int hist[hist_size] = { 0 };
			for (int i = 0; i < N; i++)
			{
				speed[i] = sqrt(particles[i].vX*particles[i].vX +
					particles[i].vY*particles[i].vY +
					particles[i].vZ*particles[i].vZ);
				hist[(int)(50 * speed[i])]++;
			}

			for (int i = 0; i < hist_size; i++)
			{
				if (max_F < hist[i])
					max_F = hist[i];
			}
			float scaleHist = 0.4 * HEIGHT / max_F;

			plotWindow.clear();
			plotWindow.draw(lineX);
			plotWindow.draw(lineY);
			for (int i = 0; i < hist_size; i++)
			{
				RectangleShape line(Vector2f(scaleHist * hist[i], 4));
				line.setPosition(21 + i * 5, HEIGHT - 100);
				line.rotate(-90);
				plotWindow.draw(line);
			}
			plotWindow.display();
		}
#endif


#ifdef DIFF
		//Diffusion part
		if (iter < diff_iters)
		{
			for (int i = 0; i < diff_part; i++)
			{
				diff_X[iter] += (particles[i].sX + particles[i].wX*Lx - init_X[i])*
					(particles[i].sX + particles[i].wX*Lx - init_X[i]);

				diff_Y[iter] += (particles[i].sY + particles[i].wY*Ly - init_Y[i])*
					(particles[i].sY + particles[i].wY*Ly - init_Y[i]);

				diff_Z[iter] += (particles[i].sZ + particles[i].wZ*Lz - init_Z[i])*
					(particles[i].sZ + particles[i].wZ*Lz - init_Z[i]);
			}
			diff_X[iter] /= diff_part;
			diff_Y[iter] /= diff_part;
			diff_Z[iter] /= diff_part;
			//cout << iter << endl;
		}
		
		if (iter == diff_iters)
		{
			/*float max_X = 0;
			float max_Y = 0;
			float max_Z = 0;
			for (int i = diff_iters - 1; i >= 0; i--)
			{
				if (max_X < diff_X[i])
					max_X = diff_X[i];
				if (max_Y < diff_Y[i])
					max_Y = diff_Y[i];
				if (max_Z < diff_Z[i])
					max_Z = diff_Z[i];
			}
			
			plotWindow.clear();
			float scaleDiffX = 0.8 * HEIGHT / max_X;
			float scaleDiffY = 0.8 * HEIGHT / max_Y;
			float scaleDiffZ = 0.8 * HEIGHT / max_Z;
			scaleDiffX = min(scaleDiffX, scaleDiffY);
			scaleDiffX = min(scaleDiffX, scaleDiffZ);
			float scaleX = 0.9 * (WIDTH - 20) / diff_iters;
			plotWindow.draw(lineX);
			lineY.setPosition(0, HEIGHT - 20);
			plotWindow.draw(lineY);
			for (int i = 0; i < diff_iters; i++)
			{
				CircleShape diffParticle(1);
				diffParticle.setFillColor(Color(255, 50, 50));
				diffParticle.setPosition(20 + i*scaleX, HEIGHT - 20 - diff_X[i] * scaleDiffX);
				plotWindow.draw(diffParticle);
				diffParticle.setFillColor(Color(50, 255, 50));
				diffParticle.setPosition(20 + i*scaleX, HEIGHT - 20 - diff_Y[i] * scaleDiffX);
				plotWindow.draw(diffParticle);
				diffParticle.setFillColor(Color(50, 50, 255));
				diffParticle.setPosition(20 + i*scaleX, HEIGHT - 20 - diff_Z[i] * scaleDiffX);
				plotWindow.draw(diffParticle);
			}
			plotWindow.display();*/

			FILE* pipe = _popen("C:/Program\" \Files/gnuplot/bin/gnuplot\" \--persist", "w");
			fprintf(pipe, "\n");
			fprintf(pipe, "load 'settings.ps'\n");
			fprintf(pipe, "set xlabel 'Time'\n");
			fprintf(pipe, "set ylabel '< x² >'\n");
			fprintf(pipe, "set key top left\n");
			fprintf(pipe, "plot '-' title 'X' with lines linestyle 1,'-' title 'Y' with lines linestyle 2,'-' title 'Z' with lines linestyle 3\n");
			for (int i = 0; i < diff_iters; i++)
				fprintf(pipe, "%f %f\n", float(i), diff_X[i]);
			fprintf(pipe, "e\n");
			for (int i = 0; i < diff_iters; i++)
				fprintf(pipe, "%f %f\n", float(i), diff_Y[i]);
			fprintf(pipe, "e\n");
			for (int i = 0; i < diff_iters; i++)
				fprintf(pipe, "%f %f\n", float(i), diff_Z[i]);
			fprintf(pipe, "e\n");
			fflush(pipe);
			fclose(pipe);

			system("pause");
			lineY.setPosition(0, HEIGHT - 100);
		}
#endif
		

#ifdef VEL_PROF
		// Velocity Distribution - Poiseuille
		if (iter < _maxIters + 50 && iter > _maxIters)
		{
			for (int i = n; i < N_tot; i++)
			{
				if (particles[i].sY < 15 && particles[i].sY > 5)
				{
					if (particles[i].sX < Lx && particles[i].sX > 0)
					{
						int index = (int)particles[i].sZ;
						velProfileX_Z[index] += particles[i].vX;
						velProfileCountX_Z[index]++;
					}
				}
			}
		}
		if (iter == _maxIters + 51)
		{
			float max_Vx = 0;
			for (int i = 0; i < Lz; i++)
			{
				if (velProfileCountX_Z[i] != 0)
					velProfileX_Z[i] /= velProfileCountX_Z[i];
				
				if (max_Vx < velProfileX_Z[i])
					max_Vx = velProfileX_Z[i];
			}
			cout << "Viscosity obtained: " << g*Lz*Lz/(8*max_Vx) << "\n\n";

			/*plotWindow.clear();
			plotWindow.draw(lineX);
			lineY.setPosition(0, HEIGHT - 20);
			plotWindow.draw(lineY);
			float scaleVelX = 0.8 * HEIGHT / max_Vx;
			float scaleX = 0.9 * (WIDTH-20) / Lz;
			
			const int nPts = 100;
			float th_VelProfileX_Z[nPts] = { 0.0 };	// 100 points for plotting theoretical velocity profile
			float deltaZ = Lz/(nPts - 1.0);
			for (int i = 0; i < nPts; i++)
			{
				float _z = i*deltaZ;
				th_VelProfileX_Z[i] = 4 * max_Vx*_z*(Lz - _z) / (Lz*Lz);
			}

			//Theoritical profile
			for (int i = 1; i < nPts; i++)
			{
				Vertex lineSegX[] =
				{
					Vertex(Vector2f(20 + (i - 1)*deltaZ*scaleX, HEIGHT - 20 - th_VelProfileX_Z[i - 1] * scaleVelX), Color(255, 150, 0)),
					Vertex(Vector2f(20 + i*deltaZ*scaleX, HEIGHT - 20 - th_VelProfileX_Z[i] * scaleVelX), Color(255, 150, 0))
				};
				plotWindow.draw(lineSegX, 2, Lines);
			}
			//Obtained profile
			for (int i = 1; i < Ly; i++)
			{
				Vertex lineSegX[] =
				{
					Vertex(Vector2f(20 + (i - 1)*scaleX, HEIGHT - 20 - velProfileX[i - 1] * scaleVelX), Color(0, 255, 0)),
					Vertex(Vector2f(20 + i * scaleX, HEIGHT - 20 - velProfileX[i] * scaleVelX), Color(0, 255, 0))
				};
				plotWindow.draw(lineSegX, 2, Lines);

				Vertex lineSegY[] =
				{
					Vertex(Vector2f(20 + (i - 1)*scaleX, HEIGHT - 20 - velProfileY[i - 1] * scaleVelX),Color(255,0,0)),
					Vertex(Vector2f(20 + i * scaleX, HEIGHT - 20 - velProfileY[i] * scaleVelX), Color(255, 0, 0))
				};
				plotWindow.draw(lineSegY, 2, Lines);
				
			}
			for (int i = 0; i < Lz; i++)
			{
				CircleShape point(2);
				point.setPosition(20 + (i+0.5)*scaleX, HEIGHT - 20 - velProfileX_Z[i] * scaleVelX);
				point.setFillColor(Color(0, 255, 0));
				plotWindow.draw(point);
			}

			plotWindow.display();*/

			FILE* pipe = _popen("C:/Program\" \Files/gnuplot/bin/gnuplot\" \--persist", "w");
			fprintf(pipe, "\n");
			fprintf(pipe, "load 'settings.ps'\n");
			fprintf(pipe, "set xlabel 'z'\n");
			fprintf(pipe, "set ylabel 'V_x (z)'\n");
			fprintf(pipe, "set key top right\n");
			fprintf(pipe, "f(x) = 4 * %f * x * (%f - x) / (%f)\n", max_Vx, float(Lz), float(Lz*Lz));
			fprintf(pipe, "plot [0:%f] f(x) title 'Theoritical' with lines linestyle 1\n",float(Lz));
			fprintf(pipe, "replot '-' title 'Obtained' with point pt 7 ps 1\n");
			for (int i = 0; i < Lz; i++)
				fprintf(pipe, "%f %f\n", float(i+0.5), velProfileX_Z[i]);
			fprintf(pipe, "e\n");
			fflush(pipe);
			fclose(pipe);

			system("pause");
			lineY.setPosition(0, HEIGHT - 100);
		}
#endif


#ifdef E_M	
		//Energy Momentum Plot
		if (iter < (EM_max_iter + EM_offset_iter) && iter >= EM_offset_iter)
		{
			energyArr[iter - EM_offset_iter] = Energy(particles, 0, n);
			momentumArr[iter - EM_offset_iter] = Momentum(particles, 0, n);
		}

		if (iter == EM_max_iter + EM_offset_iter)
		{
			/*float max_E = 0;
			float max_M = 0;
			for (int i = 0; i < EM_max_iter; i++)
			{
				if (max_E < energyArr[i])
					max_E = energyArr[i];
				if (max_M < momentumArr[i])
					max_M = momentumArr[i];
			}

			plotWindow.clear();
			plotWindow.draw(lineX);
			lineY.setPosition(0, HEIGHT - 20);
			plotWindow.draw(lineY);
			float scaleE = 0.8 * HEIGHT / max_E;
			float scaleM = 0.8 * HEIGHT / max_M;
			float scaleX = 0.9 * (WIDTH - 20) / EM_max_iter;
			for (int i = 1; i < EM_max_iter; i++)
			{
				Vertex lineE[] =
				{
					Vertex(Vector2f(20 + (i - 1)*scaleX, HEIGHT - 20 - energyArr[i - 1] * scaleE), Color(255, 0, 0)),
					Vertex(Vector2f(20 + (i)*scaleX, HEIGHT - 20 - energyArr[i] * scaleE), Color(255, 0, 0))
				};
				plotWindow.draw(lineE, 2, Lines);
				Vertex lineM[] =
				{
					Vertex(Vector2f(20 + (i - 1)*scaleX, HEIGHT - 20 - momentumArr[i - 1] * scaleM), Color(0, 255, 0)),
					Vertex(Vector2f(20 + (i)*scaleX, HEIGHT - 20 - momentumArr[i] * scaleM), Color(0, 255, 0))
				};
				plotWindow.draw(lineM, 2, Lines);
			}

			
			plotWindow.display();
			system("pause");
			lineY.setPosition(0, HEIGHT - 100);*/

			FILE *pipe = _popen("C:/Program\" \Files/gnuplot/bin/gnuplot\" \--persist", "w");
			fprintf(pipe, "\n");
			fprintf(pipe, "load 'settings.ps'\n");
			fprintf(pipe, "set xlabel 'Time'\n");
			fprintf(pipe, "set ylabel 'Energy'\n");
			fprintf(pipe, "set yrange [0:*]\n");
			fprintf(pipe, "plot '-' with lines linestyle 4\n");
			for (int i = 0; i < EM_max_iter; i++)
				fprintf(pipe, "%f %f\n", float(i), energyArr[i]);
			fflush(pipe);
			fclose(pipe);

			float avg_energy = 0.0;
			for (int i = 0; i < EM_max_iter; i++)
			{
				avg_energy += energyArr[i];
			}
			avg_energy /= EM_max_iter;
			cout << "Actual avg energy : " << avg_energy << endl;
			cout << "1.5kT : " << 1.5 * kB_T << endl;
			cout << "Energy standard deviation : " << std_dev(energyArr, avg_energy, 0, EM_max_iter - 1) << endl;


			system("pause");
		}
#endif


#if defined(E2E_LEN) && defined(POLYMER)

		if (iter < e2e_max_iter)
		{
			/*float x_ = (particles[n - 1].sX + particles[n - 1].wX*Lx - particles[0].sX - particles[0].wX*Lx);
			float y_ = (particles[n - 1].sY + particles[n - 1].wY*Ly - particles[0].sY - particles[0].wY*Ly);
			float z_ = (particles[n - 1].sZ + particles[n - 1].wZ*Lz - particles[0].sZ - particles[0].wZ*Lz);*/
			float x_ = 0, y_ = 0, z_ = 0;
			for (int i = 0; i < n; i++)
			{
				x_ += particles[i].sX + particles[i].wX*Lx;
				y_ += particles[i].sY + particles[i].wY*Ly;
				z_ += particles[i].sZ + particles[i].wZ*Lz;
			}
			float x__ = x_ / n; float y__ = y_ / n; float z__ = z_ / n; 
			for (int i = 0; i < n; i++)
			{
				x_ = particles[i].sX + particles[i].wX*Lx;
				y_ = particles[i].sY + particles[i].wY*Ly;
				z_ = particles[i].sZ + particles[i].wZ*Lz;
				e2e_len[iter] += (x_ - x__)*(x_ - x__) + (y_ - y__)*(y_ - y__) + (z_ - z__)*(z_ - z__);
			}
			e2e_len[iter] /= n;
		}
		if (iter == e2e_max_iter)
		{
			FILE *pipe = _popen("C:/Program\" \Files/gnuplot/bin/gnuplot\" \--persist", "w");
			fprintf(pipe, "\n");
			fprintf(pipe, "load 'settings.ps'\n");
			fprintf(pipe, "set xlabel 'Time'\n");
			fprintf(pipe, "set ylabel '(R_g)^2'\n");
			fprintf(pipe, "plot '-' with lines linestyle 1\n");
			for (int i = 0; i < e2e_max_iter; i++)
				fprintf(pipe, "%f %f\n", float(i), e2e_len[i]);
			fflush(pipe);
			fclose(pipe);

			system("pause");
		}

#endif
	
	}


	return 0;
}