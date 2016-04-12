// Mustafar.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define DOUBLE_PRECISION 0
#if DOUBLE_PRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

#define STEP_COUNT 10e+1
#define LAMBDA_DEN 9e+4 // 9e+5 gives epsilon out of range 

int _tmain(int argc, _TCHAR* argv[])
{
	REAL np = 92115;
	REAL G = 6.673e-11; // N m^2 / kg ^2
	REAL lambda = np / LAMBDA_DEN;
	REAL m1 = lambda * (1.9891e+30); // kg
	REAL m2 = lambda * (3.301e+23); // kg
	REAL M = m1 + m2; // kg
	REAL epsilon = pow(lambda, -1) * 0.2056; // -
	if (epsilon <= 0 || epsilon >= 1)
	{
		printf("Epsilon=%e out of range!\n", epsilon);
		getchar();
		return 0;
	}
	REAL semiMinorAxis = pow(lambda, 2) * (5.791e+10); // m
	REAL mu = G * M; // m^3 / s^2
	REAL specificAngularMomentumSquared = semiMinorAxis * mu * (1 - pow(epsilon, 2)); // m^4 / s^2

	printf("NP=%e\n", np);
	printf("LAMBDA=%e\n", lambda);
	printf("M1=%e\n", m1);
	printf("M2=%e\n", m2);
	printf("MU=%e\n", mu);
	printf("EPSILON=%e\n", epsilon);
	printf("SEMI MINOR AXIS=%e\n", semiMinorAxis);
	printf("SPECIFIC ANGULAR MOMENTUM SQUARED=%e\n", specificAngularMomentumSquared);

	REAL u_0 = pow(semiMinorAxis * (1 - epsilon), -1);
	REAL v_0 = 0;
	REAL k = (2 * M_PI) / STEP_COUNT;

	printf("U_0=%e\n", u_0);
	printf("V_0=%e\n", v_0);
	printf("STEP COUNT=%e\n", STEP_COUNT);
	printf("RADIAN STEP=%e\n", k);

	printf("Start!\n", epsilon);
	std::fstream output;
	output.open("Mustafar.csv", std::ios_base::out);

	REAL u_n = u_0;
	REAL v_n = v_0;
	REAL k_mu_sams = (k * mu * pow(specificAngularMomentumSquared, -1));
	for (long n = 0; n < STEP_COUNT; n++)
	{
		REAL theta_n = n * k;
		REAL radius_n = 1 / u_n;
		output << theta_n << ";" << radius_n << ";" << "\n";
		REAL u_n_1 = u_n + (k * v_n);
		REAL v_n_1 = v_n - (k * u_n) + k_mu_sams;
		u_n = u_n_1;
		v_n = v_n_1;
	}

	output.close();
	printf("Done!\n", epsilon);
	getchar();
	return 0;
}

