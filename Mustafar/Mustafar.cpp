// Mustafar.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

REAL computeDDR(REAL u_n, REAL semiMinorAxis, REAL epsilon)
{
	REAL q = semiMinorAxis * (1 - epsilon);
	REAL DRR = (pow(u_n, -1) / q) - 1;
	return DRR;
}

REAL algorithm1(REAL u_0, REAL v_0, REAL k, REAL alpha, long steps, std::fstream& output)
{
	REAL u_n = u_0;
	REAL v_n = v_0;
	REAL k_alpha = k * alpha;
	for (long n = 0; n < steps; n++)
	{
		REAL theta_n = n * k;
		REAL radius_n = 1 / u_n;
		output << theta_n << ";" << radius_n << ";" << "\n";

		REAL u_n_1 = u_n + (k * v_n);
		REAL v_n_1 = v_n - (k * u_n) + k_alpha;

		u_n = u_n_1;
		v_n = v_n_1;
	}
	return u_n;
}

REAL algorithm2(REAL u_0, REAL v_0, REAL k, REAL alpha, long steps, std::fstream& output)
{
	REAL u_n = u_0;
	REAL v_n = v_0;
	for (long n = 0; n < steps; n++)
	{
		REAL theta_n = n * k;
		REAL radius_n = 1 / u_n;
		output << theta_n << ";" << radius_n << ";" << "\n";

		REAL w_1 = u_n + (k * (v_n / 2));
		REAL z_1 = v_n + (k * ((alpha - u_n) / 2));
		REAL w_2 = u_n + (k * (z_1 / 2));
		REAL z_2 = v_n + (k * ((alpha - w_1) / 2));
		REAL w_3 = u_n + (k * z_2);
		REAL z_3 = v_n + (k * (alpha - w_2));

		REAL k_v_z = (k * (v_n + ((2 * z_1) + (2 * z_2) + z_3) / 6));
		REAL k_v_alfa_u_w = (k * (v_n + ((6 * alpha) - (u_n)-(2 * w_1) - (2 * w_2) - (w_3)) / 6));

		REAL u_n_1 = u_n + k_v_z;
		REAL v_n_1 = v_n + k_v_alfa_u_w;

		u_n = u_n_1;
		v_n = v_n_1;
	}
	return u_n;
}

int main(int argc, char* argv[])
{
	long long steps = 0;
	long alg = 0;
	for (int argIndex = 1; argIndex < argc; argIndex++)
	{
		if (!strcmp(argv[argIndex], "-alg") &&
			(argIndex + 1) < argc)
		{
			alg = atoll(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-steps") &&
			(argIndex + 1) < argc)
		{
			steps = atoll(argv[argIndex + 1]);
			argIndex++;
		}
	}
	
	if (steps == 0 ||
		alg <= 0
		|| alg >= 2)
	{
		printf("Wrong params!\n");
		printf("-alg \t Algorithm number {1, 2}.\n");
		printf("-steps \t Algorithm step count.\n");
		getchar();
		return 0;
	}

	char runDatabaseFileName[255];
	memset(runDatabaseFileName, 0, sizeof(runDatabaseFileName) / sizeof(runDatabaseFileName[0]));
	sprintf_s(runDatabaseFileName, "alg_%ld_steps_%lld_db.csv", alg, steps);

	printf("Creating output files...\n");
	std::fstream dbOutput;
	dbOutput.open(runDatabaseFileName, std::ios_base::out);

	std::fstream statsOutput;
	statsOutput.open("mustafar_stats.csv", std::ios_base::app);

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
	
	printf("ALGORITHM=%d\n", alg);
	statsOutput << alg << ";";
	printf("DOUBLE PRECISION=%d\n", DOUBLE_PRECISION);
	statsOutput << DOUBLE_PRECISION << ";";
	printf("NP=%e\n", np);
	statsOutput << np << ";";
	printf("LAMBDA=%e\n", lambda);
	statsOutput << lambda << ";";
	printf("M1=%e\n", m1);
	statsOutput << m1 << ";";
	printf("M2=%e\n", m2);
	statsOutput << m2 << ";";
	printf("MU=%e\n", mu);
	statsOutput << mu << ";";
	printf("EPSILON=%e\n", epsilon);
	statsOutput << epsilon << ";";
	printf("SEMI MINOR AXIS=%e\n", semiMinorAxis);
	statsOutput << semiMinorAxis << ";";
	printf("SPECIFIC ANGULAR MOMENTUM SQUARED=%e\n", specificAngularMomentumSquared);
	statsOutput << specificAngularMomentumSquared << ";";

	REAL u_0 = pow(semiMinorAxis * (1 - epsilon), -1);
	REAL v_0 = 0;
	REAL k = (2 * M_PI) / steps;

	printf("U_0=%e\n", u_0);
	statsOutput << u_0 << ";";
	printf("V_0=%e\n", v_0);
	statsOutput << v_0 << ";";
	printf("STEP COUNT=%e\n", steps);
	statsOutput << steps << ";";
	printf("RADIAN STEP=%e\n", k);
	statsOutput << k << ";";
	
	REAL alpha = mu * (pow(specificAngularMomentumSquared, -1));
	REAL u_n = 0;

	printf("Start!\n");
	auto startTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	if (alg == 1)
	{
		u_n = algorithm1(u_0, v_0, k, alpha, steps, dbOutput);
	}
	else if (alg == 2)
	{
		u_n = algorithm2(u_0, v_0, k, alpha, steps, dbOutput);
	}
	auto finishTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	auto tiempoDeCorrida = std::chrono::duration_cast<std::chrono::nanoseconds>(finishTimer - startTimer).count();
	printf("End!\n");

	printf("Saving database...\n");
	dbOutput.close();

	REAL ddr = computeDDR(u_n, semiMinorAxis, epsilon);
	printf("DDR=%e\n", ddr);

	printf("Saving stats...\n");
	statsOutput << ddr << ";" << tiempoDeCorrida << ";" << "\n";
	statsOutput.close();
	
	printf("Done!\n");
	getchar();
	return 0;
}

/*
Metodo para calcular tiempo de corrida en milisegundos si es necesario
#include <time.h>
clock_t start;
clock_t end;
double msecs;
start = clock();
end = clock();
msecs = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;
printf("Timer=%e\n", msecs);
*/