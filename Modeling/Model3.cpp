#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define PI 3.1415926535

float** space2d(int nr, int nc);
void free_space2d(float** a, int nx);
void wfile(char filename[], float** data, int nr, int nc);
void create_model(float** vp, float** vs, float** rho, float** lamda, float** mu, int nr, int nc);
float** extmodel(float** init_model, int nr, int nc, int np);

int main()
{
	int NX = 5000;		    
	int NZ = 200;		  
	int NP = 100;		  
	int sx = 2500 + NP;    
	int sz = 0;
	int NX_ext = NX + 2*NP; 
	int NZ_ext = NZ + NP;    
	int NT = 60000;		  
	int dsnap = 100;        
	double H = 1;		  
	double RC = 0.0001;	  
	double DP = NP * H;		

	double DT = 0.00005;	 
	double F0 = 30.0;		
	double T0 = 1.2 / F0;   
	double Vpmax = 6600.0;  
	double Vpmin = 3600.0;  
	double Vsmax = 2200.0;
	double Vsmin = 1200.0;

	float** vs;		
	float** vp;		
	float** rho;	
	float** mu;		
	float** lamda;
	float** sis_x;  
	float** sis_z;  

	float** vx_x;
	float** vx_z;
	float** vx;
	float** vz_x;
	float** vz_z;
	float** vz;
	float** txx_x;
	float** txx_z;
	float** tzz_x;
	float** tzz_z;
	float** txz_x;
	float** txz_z;
	float** dxi, ** dxi2;
	float** dzj, ** dzj2;

	double tt, x, z, xoleft, xoright, d0, suma;
	double v0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, rhoxz, lam2uz, mux, lamdaz, best_dt;
	int ix, iz, it;
	FILE* fp_snpx, * fp_snpz;
	float utemp, wtemp;

	vs = space2d(NZ, NX);
	vp = space2d(NZ, NX);
	rho = space2d(NZ, NX);
	mu = space2d(NZ, NX);
	lamda = space2d(NZ, NX);
	create_model(vp, vs, rho, lamda, mu, NZ, NX);

	char vs_name[] = "vs_ext.dat";
	float** vs_ext;		
	float** vp_ext;     
	float** rho_ext;	
	float** mu_ext;		
	float** lamda_ext;

	vs_ext = extmodel(vs, NZ, NX, NP);
	vp_ext = extmodel(vp, NZ, NX, NP);
	rho_ext = extmodel(rho, NZ, NX, NP);
	mu_ext = extmodel(mu, NZ, NX, NP);
	lamda_ext = extmodel(lamda, NZ, NX, NP);
	wfile(vs_name, vs_ext, NZ_ext, NX_ext);

	vx_x = space2d(NZ_ext, NX_ext);
	vx_z = space2d(NZ_ext, NX_ext);
	vx   = space2d(NZ_ext, NX_ext);
	vz_x = space2d(NZ_ext, NX_ext);
	vz_z = space2d(NZ_ext, NX_ext);
	vz   = space2d(NZ_ext, NX_ext);

	txx_x = space2d(NZ_ext, NX_ext);
	txx_z = space2d(NZ_ext, NX_ext);
	tzz_x = space2d(NZ_ext, NX_ext);
	tzz_z = space2d(NZ_ext, NX_ext);
	txz_x = space2d(NZ_ext, NX_ext);
	txz_z = space2d(NZ_ext, NX_ext);

	dxi  = space2d(NZ_ext, NX_ext);
	dxi2 = space2d(NZ_ext, NX_ext);
	dzj  = space2d(NZ_ext, NX_ext);
	dzj2 = space2d(NZ_ext, NX_ext);

	xoleft = DP;					
	xoright = (NX_ext - 1) * H - DP;	


	for (iz = 0; iz < NZ_ext; iz++)
	{
		for (ix = 0; ix < NX_ext; ix++)
		{
			x = ix * H;
			if (x < xoleft)
			{
				v0 = vp_ext[iz][0];
				d0 = 3.0 * v0 * log(1.0 / RC) / (2.0 * DP); 
				dxi[iz][ix] = d0 * pow(((xoleft - x) / DP), 2);
				dxi2[iz][ix] = d0 * pow(((xoleft - x - 0.5 * H) / DP), 2);
			}
			else if (x >= 0.9999 * xoright)
			{
				v0 = vp_ext[iz][NX_ext - 1];
				d0 = 3.0 * v0 * log(1.0 / RC) / (2.0 * DP);
				dxi[iz][ix] = d0 * pow(((x - xoright) / DP), 2);
				dxi2[iz][ix] = d0 * pow(((x + 0.5 * H - xoright) / DP), 2);
			}
			else
			{
				dxi[iz][ix] = 0.; dxi2[iz][ix] = 0.;
			}
		}
	}

	xoright = (NZ_ext - 1) * H - DP;	
	for (iz = 0; iz < NZ_ext; iz++)
	{
		z = iz * H;
		for (ix = 0; ix < NX_ext; ix++)
		{
			if (z >= 0.9999 * xoright)
			{
				v0 = vp_ext[NZ_ext - 1][ix];
				d0 = 3.0 * v0 * log(1.0 / RC) / (2.0 * DP);
				dzj[iz][ix] = d0 * pow(((z - xoright) / DP), 2);
				dzj2[iz][ix] = d0 * pow(((z + 0.5 * H - xoright) / DP), 2);
			}
			else
			{
				dzj[iz][ix] = 0.; dzj2[iz][ix] = 0.;
			}
		}
	}

	double a1 = 1.211242675781249;
	double a2 = -0.089721679687499;
	double a3 = 0.013842773437500;
	double a4 = -0.001765659877232;
	double a5 = 0.000118679470486;

	suma = fabs(a1) + fabs(a2) + fabs(a3) + fabs(a4) + fabs(a5);
	best_dt = H / (sqrt(2.0) * Vpmax * suma);
	if (DT > best_dt)
	{
		printf("Time step too large, should be smaller than %f\n", best_dt);
		exit(0);
	}

	if (Vsmin / (F0 * H) < 10)
		printf("Spatial step too large, may cause significant numerical dispersion\n");

	clock_t start, finish;
	start = clock();

	sis_x = space2d(NX, NT); sis_z = space2d(NX, NT);

	fp_snpx = fopen("vx_G.dat", "wb");
	fp_snpz = fopen("vz_G.dat", "wb");
	fclose(fp_snpx);
	fclose(fp_snpz);

	for (it = 0; it < NT; it++)
	{
		tt = it * DT;
		if (it % 100 == 0) printf("it=%d\n", it);

		for (iz = 0; iz < NZ_ext - 5; iz++)
			for (ix = 5; ix < NX_ext - 5; ix++)
			{
				rhoxz = 1.0 / (0.25 * (1 / rho_ext[iz][ix] + 1 / rho_ext[iz][ix + 1] + 1 / rho_ext[iz + 1][ix] + 1 / rho_ext[iz + 1][ix + 1]));
				y1 = (DT / (rhoxz * H)) * (a1 * (txx_x[iz][ix + 1] - txx_x[iz][ix] + txx_z[iz][ix + 1] - txx_z[iz][ix]) +
					a2 * (txx_x[iz][ix + 2] - txx_x[iz][ix - 1] + txx_z[iz][ix + 2] - txx_z[iz][ix - 1]) +
					a3 * (txx_x[iz][ix + 3] - txx_x[iz][ix - 2] + txx_z[iz][ix + 3] - txx_z[iz][ix - 2]) +
					a4 * (txx_x[iz][ix + 4] - txx_x[iz][ix - 3] + txx_z[iz][ix + 4] - txx_z[iz][ix - 3]) +
					a5 * (txx_x[iz][ix + 5] - txx_x[iz][ix - 4] + txx_z[iz][ix + 5] - txx_z[iz][ix - 4]));
				if (iz == 0)
					y2 = (DT / (rhoxz * H)) * (a1 * (txz_x[iz + 1][ix] - txz_x[iz][ix] + txz_z[iz + 1][ix] - txz_z[iz][ix]) +
						a2 * (txz_x[iz + 2][ix] + txz_x[iz + 1][ix] + txz_z[iz + 2][ix] + txz_z[iz + 1][ix]) +
						a3 * (txz_x[iz + 3][ix] + txz_x[iz + 2][ix] + txz_z[iz + 3][ix] + txz_z[iz + 2][ix]) +
						a4 * (txz_x[iz + 4][ix] + txz_x[iz + 3][ix] + txz_z[iz + 4][ix] + txz_z[iz + 3][ix]) +
						a5 * (txz_x[iz + 5][ix] + txz_x[iz + 4][ix] + txz_z[iz + 5][ix] + txz_z[iz + 4][ix]));
				else if (iz == 1)
					y2 = (DT / (rhoxz * H)) * (a1 * (txz_x[iz + 1][ix] - txz_x[iz][ix] + txz_z[iz + 1][ix] - txz_z[iz][ix]) +
						a2 * (txz_x[iz + 2][ix] - txz_x[iz - 1][ix] + txz_z[iz + 2][ix] - txz_z[iz - 1][ix]) +
						a3 * (txz_x[iz + 3][ix] + txz_x[iz + 2][ix] + txz_z[iz + 3][ix] + txz_z[iz + 2][ix]) +
						a4 * (txz_x[iz + 4][ix] + txz_x[iz + 3][ix] + txz_z[iz + 4][ix] + txz_z[iz + 3][ix]) +
						a5 * (txz_x[iz + 5][ix] + txz_x[iz + 4][ix] + txz_z[iz + 5][ix] + txz_z[iz + 4][ix]));
				else if (iz == 2)
					y2 = (DT / (rhoxz * H)) * (a1 * (txz_x[iz + 1][ix] - txz_x[iz][ix] + txz_z[iz + 1][ix] - txz_z[iz][ix]) +
						a2 * (txz_x[iz + 2][ix] - txz_x[iz - 1][ix] + txz_z[iz + 2][ix] - txz_z[iz - 1][ix]) +
						a3 * (txz_x[iz + 3][ix] - txz_x[iz - 2][ix] + txz_z[iz + 3][ix] - txz_z[iz - 2][ix]) +
						a4 * (txz_x[iz + 4][ix] + txz_x[iz + 3][ix] + txz_z[iz + 4][ix] + txz_z[iz + 3][ix]) +
						a5 * (txz_x[iz + 5][ix] + txz_x[iz + 4][ix] + txz_z[iz + 5][ix] + txz_z[iz + 4][ix]));
				else if (iz == 3)
					y2 = (DT / (rhoxz * H)) * (a1 * (txz_x[iz + 1][ix] - txz_x[iz][ix] + txz_z[iz + 1][ix] - txz_z[iz][ix]) +
						a2 * (txz_x[iz + 2][ix] - txz_x[iz - 1][ix] + txz_z[iz + 2][ix] - txz_z[iz - 1][ix]) +
						a3 * (txz_x[iz + 3][ix] - txz_x[iz - 2][ix] + txz_z[iz + 3][ix] - txz_z[iz - 2][ix]) +
						a4 * (txz_x[iz + 4][ix] - txz_x[iz - 3][ix] + txz_z[iz + 4][ix] - txz_z[iz - 3][ix]) +
						a5 * (txz_x[iz + 5][ix] + txz_x[iz + 4][ix] + txz_z[iz + 5][ix] + txz_z[iz + 4][ix]));
				else
					y2 = (DT / (rhoxz * H)) * (a1 * (txz_x[iz + 1][ix] - txz_x[iz][ix] + txz_z[iz + 1][ix] - txz_z[iz][ix]) +
						a2 * (txz_x[iz + 2][ix] - txz_x[iz - 1][ix] + txz_z[iz + 2][ix] - txz_z[iz - 1][ix]) +
						a3 * (txz_x[iz + 3][ix] - txz_x[iz - 2][ix] + txz_z[iz + 3][ix] - txz_z[iz - 2][ix]) +
						a4 * (txz_x[iz + 4][ix] - txz_x[iz - 3][ix] + txz_z[iz + 4][ix] - txz_z[iz - 3][ix]) +
						a5 * (txz_x[iz + 5][ix] - txz_x[iz - 4][ix] + txz_z[iz + 5][ix] - txz_z[iz - 4][ix]));
				y3 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (txz_x[iz][ix] - txz_x[iz][ix - 1] + txz_z[iz][ix] - txz_z[iz][ix - 1]) +
					a2 * (txz_x[iz][ix + 1] - txz_x[iz][ix - 2] + txz_z[iz][ix + 1] - txz_z[iz][ix - 2]) +
					a3 * (txz_x[iz][ix + 2] - txz_x[iz][ix - 3] + txz_z[iz][ix + 2] - txz_z[iz][ix - 3]) +
					a4 * (txz_x[iz][ix + 3] - txz_x[iz][ix - 4] + txz_z[iz][ix + 3] - txz_z[iz][ix - 4]) +
					a5 * (txz_x[iz][ix + 4] - txz_x[iz][ix - 5] + txz_z[iz][ix + 4] - txz_z[iz][ix - 5]));
				if (iz == 0)
					y4 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (tzz_x[iz][ix] + tzz_x[iz + 1][ix] + tzz_z[iz][ix] + tzz_z[iz + 1][ix]) +
						a2 * (tzz_x[iz + 1][ix] + tzz_x[iz + 2][ix] + tzz_z[iz + 1][ix] + tzz_z[iz + 2][ix]) +
						a3 * (tzz_x[iz + 2][ix] + tzz_x[iz + 3][ix] + tzz_z[iz + 2][ix] + tzz_z[iz + 3][ix]) +
						a4 * (tzz_x[iz + 3][ix] + tzz_x[iz + 4][ix] + tzz_z[iz + 3][ix] + tzz_z[iz + 4][ix]) +
						a5 * (tzz_x[iz + 4][ix] + tzz_x[iz + 5][ix] + tzz_z[iz + 4][ix] + tzz_z[iz + 5][ix]));
				else if (iz == 1)
					y4 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (tzz_x[iz][ix] - tzz_x[iz - 1][ix] + tzz_z[iz][ix] - tzz_z[iz - 1][ix]) +
						a2 * (tzz_x[iz + 1][ix] + tzz_x[iz + 2][ix] + tzz_z[iz + 1][ix] + tzz_z[iz + 2][ix]) +
						a3 * (tzz_x[iz + 2][ix] + tzz_x[iz + 3][ix] + tzz_z[iz + 2][ix] + tzz_z[iz + 3][ix]) +
						a4 * (tzz_x[iz + 3][ix] + tzz_x[iz + 4][ix] + tzz_z[iz + 3][ix] + tzz_z[iz + 4][ix]) +
						a5 * (tzz_x[iz + 4][ix] + tzz_x[iz + 5][ix] + tzz_z[iz + 4][ix] + tzz_z[iz + 5][ix]));
				else if (iz == 2)
					y4 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (tzz_x[iz][ix] - tzz_x[iz - 1][ix] + tzz_z[iz][ix] - tzz_z[iz - 1][ix]) +
						a2 * (tzz_x[iz + 1][ix] - tzz_x[iz - 2][ix] + tzz_z[iz + 1][ix] - tzz_z[iz - 2][ix]) +
						a3 * (tzz_x[iz + 2][ix] + tzz_x[iz + 3][ix] + tzz_z[iz + 2][ix] + tzz_z[iz + 3][ix]) +
						a4 * (tzz_x[iz + 3][ix] + tzz_x[iz + 4][ix] + tzz_z[iz + 3][ix] + tzz_z[iz + 4][ix]) +
						a5 * (tzz_x[iz + 4][ix] + tzz_x[iz + 5][ix] + tzz_z[iz + 4][ix] + tzz_z[iz + 5][ix]));
				else if (iz == 3)
					y4 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (tzz_x[iz][ix] - tzz_x[iz - 1][ix] + tzz_z[iz][ix] - tzz_z[iz - 1][ix]) +
						a2 * (tzz_x[iz + 1][ix] - tzz_x[iz - 2][ix] + tzz_z[iz + 1][ix] - tzz_z[iz - 2][ix]) +
						a3 * (tzz_x[iz + 2][ix] - tzz_x[iz - 3][ix] + tzz_z[iz + 2][ix] - tzz_z[iz - 3][ix]) +
						a4 * (tzz_x[iz + 3][ix] + tzz_x[iz + 4][ix] + tzz_z[iz + 3][ix] + tzz_z[iz + 4][ix]) +
						a5 * (tzz_x[iz + 4][ix] + tzz_x[iz + 5][ix] + tzz_z[iz + 4][ix] + tzz_z[iz + 5][ix]));
				else if (iz == 4)
					y4 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (tzz_x[iz][ix] - tzz_x[iz - 1][ix] + tzz_z[iz][ix] - tzz_z[iz - 1][ix]) +
						a2 * (tzz_x[iz + 1][ix] - tzz_x[iz - 2][ix] + tzz_z[iz + 1][ix] - tzz_z[iz - 2][ix]) +
						a3 * (tzz_x[iz + 2][ix] - tzz_x[iz - 3][ix] + tzz_z[iz + 2][ix] - tzz_z[iz - 3][ix]) +
						a4 * (tzz_x[iz + 3][ix] - tzz_x[iz - 4][ix] + tzz_z[iz + 3][ix] - tzz_z[iz - 4][ix]) +
						a5 * (tzz_x[iz + 4][ix] + tzz_x[iz + 5][ix] + tzz_z[iz + 4][ix] + tzz_z[iz + 5][ix]));

				else
					y4 = (DT / (rho_ext[iz][ix] * H)) * (a1 * (tzz_x[iz][ix] - tzz_x[iz - 1][ix] + tzz_z[iz][ix] - tzz_z[iz - 1][ix]) +
						a2 * (tzz_x[iz + 1][ix] - tzz_x[iz - 2][ix] + tzz_z[iz + 1][ix] - tzz_z[iz - 2][ix]) +
						a3 * (tzz_x[iz + 2][ix] - tzz_x[iz - 3][ix] + tzz_z[iz + 2][ix] - tzz_z[iz - 3][ix]) +
						a4 * (tzz_x[iz + 3][ix] - tzz_x[iz - 4][ix] + tzz_z[iz + 3][ix] - tzz_z[iz - 4][ix]) +
						a5 * (tzz_x[iz + 4][ix] - tzz_x[iz - 5][ix] + tzz_z[iz + 4][ix] - tzz_z[iz - 5][ix]));

				vx_x[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dxi2[iz][ix])) * ((1 - 0.5 * DT * dxi2[iz][ix]) * vx_x[iz][ix] + y1);
				vx_z[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dzj2[iz][ix])) * ((1 - 0.5 * DT * dzj2[iz][ix]) * vx_z[iz][ix] + y2);
				vz_x[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dxi[iz][ix])) * ((1 - 0.5 * DT * dxi[iz][ix]) * vz_x[iz][ix] + y3);
				vz_z[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dzj[iz][ix])) * ((1 - 0.5 * DT * dzj[iz][ix]) * vz_z[iz][ix] + y4);
			}

		for (iz = 0; iz < NZ_ext - 5; iz++)
			for (ix = 5; ix < NX_ext - 5; ix++)
			{
				lam2uz = 0.5 * ((lamda_ext[iz][ix] + 2 * mu_ext[iz][ix]) + (lamda_ext[iz + 1][ix] + 2 * mu_ext[iz + 1][ix]));
				lamdaz = 0.5 * (lamda_ext[iz][ix] + lamda_ext[iz + 1][ix]);
				mux = 0.5 * (mu_ext[iz][ix] + mu_ext[iz][ix + 1]);

				y5 = (DT * lam2uz / H) * (a1 * (vx_x[iz][ix] - vx_x[iz][ix - 1] + vx_z[iz][ix] - vx_z[iz][ix - 1]) +
					a2 * (vx_x[iz][ix + 1] - vx_x[iz][ix - 2] + vx_z[iz][ix + 1] - vx_z[iz][ix - 2]) +
					a3 * (vx_x[iz][ix + 2] - vx_x[iz][ix - 3] + vx_z[iz][ix + 2] - vx_z[iz][ix - 3]) +
					a4 * (vx_x[iz][ix + 3] - vx_x[iz][ix - 4] + vx_z[iz][ix + 3] - vx_z[iz][ix - 4]) +
					a5 * (vx_x[iz][ix + 4] - vx_x[iz][ix - 5] + vx_z[iz][ix + 4] - vx_z[iz][ix - 5]));
				if (iz == 0)
					y6 = (lamdaz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] + vz_z[iz + 2][ix]) +
						a3 * (vz_x[iz + 3][ix] + vz_z[iz + 3][ix]) +
						a4 * (vz_x[iz + 4][ix] + vz_z[iz + 4][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else if (iz == 1)
					y6 = (lamdaz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] + vz_z[iz + 3][ix]) +
						a4 * (vz_x[iz + 4][ix] + vz_z[iz + 4][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else if (iz == 2)
					y6 = (lamdaz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] - vz_x[iz - 2][ix] + vz_z[iz + 3][ix] - vz_z[iz - 2][ix]) +
						a4 * (vz_x[iz + 4][ix] + vz_z[iz + 4][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else if (iz == 3)
					y6 = (lamdaz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] - vz_x[iz - 2][ix] + vz_z[iz + 3][ix] - vz_z[iz - 2][ix]) +
						a4 * (vz_x[iz + 4][ix] - vz_x[iz - 3][ix] + vz_z[iz + 4][ix] - vz_z[iz - 3][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else
					y6 = (lamdaz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] - vz_x[iz - 2][ix] + vz_z[iz + 3][ix] - vz_z[iz - 2][ix]) +
						a4 * (vz_x[iz + 4][ix] - vz_x[iz - 3][ix] + vz_z[iz + 4][ix] - vz_z[iz - 3][ix]) +
						a5 * (vz_x[iz + 5][ix] - vz_x[iz - 4][ix] + vz_z[iz + 5][ix] - vz_z[iz - 4][ix]));
				y7 = (lamdaz * DT / H) * (a1 * (vx_x[iz][ix] - vx_x[iz][ix - 1] + vx_z[iz][ix] - vx_z[iz][ix - 1]) +
					a2 * (vx_x[iz][ix + 1] - vx_x[iz][ix - 2] + vx_z[iz][ix + 1] - vx_z[iz][ix - 2]) +
					a3 * (vx_x[iz][ix + 2] - vx_x[iz][ix - 3] + vx_z[iz][ix + 2] - vx_z[iz][ix - 3]) +
					a4 * (vx_x[iz][ix + 3] - vx_x[iz][ix - 4] + vx_z[iz][ix + 3] - vx_z[iz][ix - 4]) +
					a5 * (vx_x[iz][ix + 4] - vx_x[iz][ix - 5] + vx_z[iz][ix + 4] - vx_z[iz][ix - 5]));
				if (iz == 0)
					y8 = (lam2uz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] + vz_z[iz + 2][ix]) +
						a3 * (vz_x[iz + 3][ix] + vz_z[iz + 3][ix]) +
						a4 * (vz_x[iz + 4][ix] + vz_z[iz + 4][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else if (iz == 1)
					y8 = (lam2uz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] + vz_z[iz + 3][ix]) +
						a4 * (vz_x[iz + 4][ix] + vz_z[iz + 4][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else if (iz == 2)
					y8 = (lam2uz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] - vz_x[iz - 2][ix] + vz_z[iz + 3][ix] - vz_z[iz - 2][ix]) +
						a4 * (vz_x[iz + 4][ix] + vz_z[iz + 4][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else if (iz == 3)
					y8 = (lam2uz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] - vz_x[iz - 2][ix] + vz_z[iz + 3][ix] - vz_z[iz - 2][ix]) +
						a4 * (vz_x[iz + 4][ix] - vz_x[iz - 3][ix] + vz_z[iz + 4][ix] - vz_z[iz - 3][ix]) +
						a5 * (vz_x[iz + 5][ix] + vz_z[iz + 5][ix]));
				else
					y8 = (lam2uz * DT / H) * (a1 * (vz_x[iz + 1][ix] - vz_x[iz][ix] + vz_z[iz + 1][ix] - vz_z[iz][ix]) +
						a2 * (vz_x[iz + 2][ix] - vz_x[iz - 1][ix] + vz_z[iz + 2][ix] - vz_z[iz - 1][ix]) +
						a3 * (vz_x[iz + 3][ix] - vz_x[iz - 2][ix] + vz_z[iz + 3][ix] - vz_z[iz - 2][ix]) +
						a4 * (vz_x[iz + 4][ix] - vz_x[iz - 3][ix] + vz_z[iz + 4][ix] - vz_z[iz - 3][ix]) +
						a5 * (vz_x[iz + 5][ix] - vz_x[iz - 4][ix] + vz_z[iz + 5][ix] - vz_z[iz - 4][ix]));
				y9 = (mux * DT / H) * (a1 * (vz_x[iz][ix + 1] - vz_x[iz][ix] + vz_z[iz][ix + 1] - vz_z[iz][ix]) +
					a2 * (vz_x[iz][ix + 2] - vz_x[iz][ix - 1] + vz_z[iz][ix + 2] - vz_z[iz][ix - 1]) +
					a3 * (vz_x[iz][ix + 3] - vz_x[iz][ix - 2] + vz_z[iz][ix + 3] - vz_z[iz][ix - 2]) +
					a4 * (vz_x[iz][ix + 4] - vz_x[iz][ix - 3] + vz_z[iz][ix + 4] - vz_z[iz][ix - 3]) +
					a5 * (vz_x[iz][ix + 5] - vz_x[iz][ix - 4] + vz_z[iz][ix + 5] - vz_z[iz][ix - 4]));
				if (iz == 0)
					y10 = (mux * DT / H) * (a1 * (vx_x[iz][ix] + vx_z[iz][ix]) +
						a2 * (vx_x[iz + 1][ix] + vx_z[iz + 1][ix]) +
						a3 * (vx_x[iz + 2][ix] + vx_z[iz + 2][ix]) +
						a4 * (vx_x[iz + 3][ix] + vx_z[iz + 3][ix]) +
						a5 * (vx_x[iz + 4][ix] + vx_z[iz + 4][ix]));
				else if (iz == 1)
					y10 = (mux * DT / H) * (a1 * (vx_x[iz][ix] - vx_x[iz - 1][ix] + vx_z[iz][ix] - vx_z[iz - 1][ix]) +
						a2 * (vx_x[iz + 1][ix] + vx_z[iz + 1][ix]) +
						a3 * (vx_x[iz + 2][ix] + vx_z[iz + 2][ix]) +
						a4 * (vx_x[iz + 3][ix] + vx_z[iz + 3][ix]) +
						a5 * (vx_x[iz + 4][ix] + vx_z[iz + 4][ix]));
				else if (iz == 2)
					y10 = (mux * DT / H) * (a1 * (vx_x[iz][ix] - vx_x[iz - 1][ix] + vx_z[iz][ix] - vx_z[iz - 1][ix]) +
						a2 * (vx_x[iz + 1][ix] - vx_x[iz - 2][ix] + vx_z[iz + 1][ix] - vx_z[iz - 2][ix]) +
						a3 * (vx_x[iz + 2][ix] + vx_z[iz + 2][ix]) +
						a4 * (vx_x[iz + 3][ix] + vx_z[iz + 3][ix]) +
						a5 * (vx_x[iz + 4][ix] + vx_z[iz + 4][ix]));
				else if (iz == 3)
					y10 = (mux * DT / H) * (a1 * (vx_x[iz][ix] - vx_x[iz - 1][ix] + vx_z[iz][ix] - vx_z[iz - 1][ix]) +
						a2 * (vx_x[iz + 1][ix] - vx_x[iz - 2][ix] + vx_z[iz + 1][ix] - vx_z[iz - 2][ix]) +
						a3 * (vx_x[iz + 2][ix] - vx_x[iz - 3][ix] + vx_z[iz + 2][ix] - vx_z[iz - 3][ix]) +
						a4 * (vx_x[iz + 3][ix] + vx_z[iz + 3][ix]) +
						a5 * (vx_x[iz + 4][ix] + vx_z[iz + 4][ix]));
				else if (iz == 4)
					y10 = (mux * DT / H) * (a1 * (vx_x[iz][ix] - vx_x[iz - 1][ix] + vx_z[iz][ix] - vx_z[iz - 1][ix]) +
						a2 * (vx_x[iz + 1][ix] - vx_x[iz - 2][ix] + vx_z[iz + 1][ix] - vx_z[iz - 2][ix]) +
						a3 * (vx_x[iz + 2][ix] - vx_x[iz - 3][ix] + vx_z[iz + 2][ix] - vx_z[iz - 3][ix]) +
						a4 * (vx_x[iz + 3][ix] - vx_x[iz - 4][ix] + vx_z[iz + 3][ix] - vx_z[iz - 4][ix]) +
						a5 * (vx_x[iz + 4][ix] + vx_z[iz + 4][ix]));
				else
					y10 = (mux * DT / H) * (a1 * (vx_x[iz][ix] - vx_x[iz - 1][ix] + vx_z[iz][ix] - vx_z[iz - 1][ix]) +
						a2 * (vx_x[iz + 1][ix] - vx_x[iz - 2][ix] + vx_z[iz + 1][ix] - vx_z[iz - 2][ix]) +
						a3 * (vx_x[iz + 2][ix] - vx_x[iz - 3][ix] + vx_z[iz + 2][ix] - vx_z[iz - 3][ix]) +
						a4 * (vx_x[iz + 3][ix] - vx_x[iz - 4][ix] + vx_z[iz + 3][ix] - vx_z[iz - 4][ix]) +
						a5 * (vx_x[iz + 4][ix] - vx_x[iz - 5][ix] + vx_z[iz + 4][ix] - vx_z[iz - 5][ix]));

				txx_x[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dxi[iz][ix])) * ((1 - 0.5 * DT * dxi[iz][ix]) * txx_x[iz][ix] + y5);
				txx_z[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dzj2[iz][ix])) * ((1 - 0.5 * DT * dzj2[iz][ix]) * txx_z[iz][ix] + y6);
				tzz_x[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dxi[iz][ix])) * ((1 - 0.5 * DT * dxi[iz][ix]) * tzz_x[iz][ix] + y7);
				tzz_z[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dzj2[iz][ix])) * ((1 - 0.5 * DT * dzj2[iz][ix]) * tzz_z[iz][ix] + y8);

				if (iz == 0)
				{
					txz_x[iz][ix] = 0; txz_z[iz][ix] = 0;
				}
				else
				{
					txz_x[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dxi2[iz][ix])) * ((1 - 0.5 * DT * dxi2[iz][ix]) * txz_x[iz][ix] + y9);
					txz_z[iz][ix] = (1.0 / (1.0 + 0.5 * DT * dzj[iz][ix])) * ((1 - 0.5 * DT * dzj[iz][ix]) * txz_z[iz][ix] + y10);
				}
			}

		vz_x[sz][sx] = vz_x[sz][sx] - 2 * PI * PI * F0 * F0 * (tt - T0) * exp(-PI * PI * F0 * F0 * (tt - T0) * (tt - T0));

		for (iz = 0; iz < NZ_ext; iz++)
		{
			vx_x[iz][0] = 0.; vx_x[iz][NX_ext - 1] = 0.; vx_z[iz][0] = 0.; vx_z[iz][NX_ext - 1] = 0.;
			vz_x[iz][0] = 0.; vz_x[iz][NX_ext - 1] = 0.; vz_z[iz][0] = 0.; vz_z[iz][NX_ext - 1] = 0.;
		}
		for (ix = 0; ix < NX_ext; ix++)
		{
			vx_x[NZ_ext - 1][ix] = 0.; vx_z[NZ_ext - 1][ix] = 0.;
			vz_x[NZ_ext - 1][ix] = 0.; vz_z[NZ_ext - 1][ix] = 0.;
		}

		for (ix = NP; ix < NP + NX; ix++)
		{
			sis_z[ix - NP][it] = vz_x[0][ix] + vz_z[0][ix];
			sis_x[ix - NP][it] = vx_x[0][ix] + vx_z[0][ix];
		}
	}

	finish = clock();
	printf("%f seconds\n", (double)(finish - start) / CLOCKS_PER_SEC);

	char sisxname[] = "sisx_G.dat";
	char siszname[] = "sisz_G.dat";
	wfile(sisxname, sis_x, NX, NT);
	wfile(siszname, sis_z, NX, NT);

	free_space2d(tzz_x, NZ_ext);
	free_space2d(tzz_z, NZ_ext);
	free_space2d(txx_x, NZ_ext);
	free_space2d(txx_z, NZ_ext);
	free_space2d(txz_x, NZ_ext);
	free_space2d(txz_z, NZ_ext);
	free_space2d(vx_x, NZ_ext);
	free_space2d(vx_z, NZ_ext);
	free_space2d(vx, NZ_ext);
	free_space2d(vz_x, NZ_ext);
	free_space2d(vz_z, NZ_ext);
	free_space2d(vz, NZ_ext);
	free_space2d(vs_ext, NZ_ext);
	free_space2d(vp_ext, NZ_ext);
	free_space2d(rho_ext, NZ_ext);
	free_space2d(vs, NZ);
	free_space2d(vp, NZ);
	free_space2d(mu, NZ);
	free_space2d(mu_ext, NZ_ext);
	free_space2d(lamda, NZ);
	free_space2d(lamda_ext, NZ_ext);
	free_space2d(rho, NZ);
	free_space2d(sis_x, NX);
	free_space2d(sis_z, NX);
	free_space2d(dxi, NZ_ext);
	free_space2d(dxi2, NZ_ext);
	free_space2d(dzj, NZ_ext);
	free_space2d(dzj2, NZ_ext);

	return 0;
}

float** space2d(int nr, int nc)
{
	float** a;
	int i;
	a = (float**)calloc(nr, sizeof(float*));
	for (i = 0; i < nr; i++)
		a[i] = (float*)calloc(nc, sizeof(float));
	return a;
}

void free_space2d(float** a, int nr)
{
	int i;
	for (i = 0; i < nr; i++)
		free(a[i]);
	free(a);
}

void wfile(char filename[], float** data, int nr, int nc)
{
	int i, j;
	FILE* fp = fopen(filename, "wb");
	for (i = 0; i < nr; i++)
		for (j = 0; j < nc; j++)
			fwrite(&data[i][j], 1, sizeof(float), fp);
	fclose(fp);
}

void create_model(float** vp, float** vs, float** rho, float** lamda, float** mu, int nr, int nc)
{
	int ix, iz;
	double lam2u;
	for (iz = 0; iz < nr; iz++)
		for (ix = 0; ix < nc; ix++)
		{
			if (iz <= 15)
			{
				vp[iz][ix] = 3600.0;
				vs[iz][ix] = 1200.0;
				rho[iz][ix] = 2000.0;
			}
			else if (iz > 15 && iz <= 50)
			{
				vp[iz][ix] = 5400.0;
				vs[iz][ix] = 1800.0;
				rho[iz][ix] = 2000.0;
			}
			else if (iz > 50 && iz <= 80)
			{
				vp[iz][ix] = 4200.0;
				vs[iz][ix] = 1400.0;
				rho[iz][ix] = 2000.0;
			}
			else
			{
				vp[iz][ix] = 6600.0;
				vs[iz][ix] = 2200.0;
				rho[iz][ix] = 2000.0;
			}
			mu[iz][ix] = rho[iz][ix] * vs[iz][ix] * vs[iz][ix];
			lam2u = rho[iz][ix] * vp[iz][ix] * vp[iz][ix];
			lamda[iz][ix] = lam2u - 2.0 * mu[iz][ix];
		}
}

float** extmodel(float** init_model, int nz, int nx, int np)
{
	float** p;
	int i, j;
	int nx2 = nx + 2 * np;
	int nz2 = nz + np;

	p = space2d(nz2, nx2);

	for (i = 0; i < nz; i++)
		for (j = 0; j < np; j++)
			p[i][j] = init_model[i][0];
	for (i = 0; i < nz; i++)
		for (j = nx; j < nx2; j++)
			p[i][j] = init_model[i][nx - 1];
	for (i = nz; i < nz2; i++)
		for (j = np; j < np + nx; j++)
			p[i][j] = init_model[nz - 1][j];
	for (i = nz; i < nz2; i++)
		for (j = 0; j < np; j++)
			p[i][j] = init_model[nz - 1][0];
	for (i = nz; i < nz2; i++)
		for (j = nx; j < nx2; j++)
			p[i][j] = init_model[nz - 1][nx - 1];
	for (i = 0; i < nz; i++)
		for (j = np; j < nx + np; j++)
			p[i][j] = init_model[i][j - np];

	return p;
}
