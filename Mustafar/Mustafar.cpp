/*
**	Universidad de Buenos Aires
**	Facultad de Ingeniería
**	75.12 Análisis Numérico I
**	Trabajo Práctico I
**
**	Merlo Leiva Nahuel –   Fabrizio Luis Cozza
**	Padrón 92115       -   Padrón 97402
*/

#include "stdafx.h"
#include "D2D1Graph.h"

template<typename T>
T algorithm1GR(T u_0, T v_0, T k, T alpha, T beta, double steps, CD2D1Graph* pGraph = nullptr)
{
	T u_n = u_0;
	T v_n = v_0;
	T k_alpha = k * alpha;
	T k_beta = k * beta;
	for (double n = 0; n < steps; n++)
	{
		if (pGraph)
		{
			REAL theta_n = n * k;
			REAL radius_n = 1 / u_n;
			pGraph->DrawPointPolar(radius_n, theta_n);
		}

		T u_n_1 = u_n + (k * v_n);
		T v_n_1 = v_n - (k * u_n) + k_alpha + k_beta * (pow(u_n, 2)); // u_n^2

		u_n = u_n_1;
		v_n = v_n_1;
	}

	if (pGraph)
	{
		REAL theta_n = steps * k;
		REAL radius_n = 1 / u_n;
		pGraph->DrawPointPolar(radius_n, theta_n);
	}
	return u_n;
}

template<typename T>
T algorithm2(T u_0, T v_0, T k, T alpha, double steps, CD2D1Graph* pGraph = nullptr)
{
	T u_n = u_0;
	T v_n = v_0;
	for (double n = 0; n < steps; n++)
	{
		if (pGraph)
		{
			REAL theta_n = n * k;
			REAL radius_n = 1 / u_n;
			pGraph->DrawPointPolar(radius_n, theta_n);
		}

		T w_1 = u_n + (k * v_n) / 2;
		T z_1 = v_n + (k * (alpha - u_n)) / 2;
		T w_2 = u_n + (k * z_1) / 2;
		T z_2 = v_n + (k * (alpha - w_1)) / 2;
		T w_3 = u_n + (k * z_2);
		T z_3 = v_n + (k * (alpha - w_2));

		T u_n_1 = u_n + (k * (v_n + ((2 * z_1) + (2 * z_2) + z_3))) / 6;
		T v_n_1 = v_n + (k * ((6 * alpha) - u_n - (2 * w_1) - (2 * w_2) - w_3)) / 6;

		u_n = u_n_1;
		v_n = v_n_1;
	}

	if (pGraph)
	{
		REAL theta_n = steps * k;
		REAL radius_n = 1 / u_n;
		pGraph->DrawPointPolar(radius_n, theta_n);
	}
	return u_n;
}

int main(int argc, char* argv[])
{
	double steps = 0;
	long alg = 0;
	long dp = 0;
	long print = 0;
	double scale = 0;
	double width = 0;
	double height = 0;
	for (int argIndex = 1; argIndex < argc; argIndex++)
	{
		if (!strcmp(argv[argIndex], "-alg") &&
			(argIndex + 1) < argc)
		{
			alg = atol(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-steps") &&
			(argIndex + 1) < argc)
		{
			steps = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-dp") &&
			(argIndex + 1) < argc)
		{
			dp = atol(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-print") &&
			(argIndex + 1) < argc)
		{
			print = atol(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-width") &&
			(argIndex + 1) < argc)
		{
			width = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-height") &&
			(argIndex + 1) < argc)
		{
			height = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-scale") &&
			(argIndex + 1) < argc)
		{
			scale = atof(argv[argIndex + 1]);
			argIndex++;
		}
	}
	
	if (steps == 0 ||
		alg <= 0
		|| alg > 2)
	{
		printf("Wrong params!\n");
		printf("-alg \t Algorithm number {1, 2}.\n");
		printf("-steps \t Algorithm step count.\n");
		printf("-dp \t Use double precision. {0, 1}\n");
		printf("-print \t Print results. {0, 1}\n");
		printf("-scale \t Print scale.\n");
		printf("-width \t Print width.\n");
		printf("-height \t Print height.\n");
		return 0;
	}

	char runDatabaseFileName[255];
	memset(runDatabaseFileName, 0, sizeof(runDatabaseFileName) / sizeof(runDatabaseFileName[0]));
	sprintf_s(runDatabaseFileName, "alg_%ld_steps_%.0f_dp_%ld_db.csv", alg, steps, dp);

	printf("Creating output files...\n");

	CD2D1Graph* pGraph = nullptr;
	if (print)
	{
		pGraph = new CD2D1Graph(width, height, scale);
		pGraph->BeginDraw();
	}

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
		return 0;
	}
	REAL semiMinorAxis = pow(lambda, 2) * (5.791e+10); // m
	REAL mu = G * M; // m^3 / s^2
	REAL specificAngularMomentumSquared = semiMinorAxis * mu * (1 - pow(epsilon, 2)); // m^4 / s^2
	REAL u_0 = pow(semiMinorAxis * (1 - epsilon), -1);
	REAL v_0 = 0;
	REAL k = (2 * M_PI) / steps;
	REAL lightningSpeed = 3 * 10e8;
	REAL lightningSpeedSquared = pow(lightningSpeed, 2); // c^2

	printf("ALGORITHM=%d\n", alg);
	printf("DOUBLE PRECISION=%d\n", dp);
	printf("NP=%e\n", np);
	printf("LAMBDA=%e\n", lambda);
	printf("M1=%e\n", m1);
	printf("M2=%e\n", m2);
	printf("MU=%e\n", mu);
	printf("EPSILON=%e\n", epsilon);
	printf("SEMI MINOR AXIS=%e\n", semiMinorAxis);
	printf("SPECIFIC ANGULAR MOMENTUM SQUARED=%e\n", specificAngularMomentumSquared);
	printf("LIGHTNING SPEED=%e\n", lightningSpeedSquared);
	printf("U_0=%e\n", u_0);
	printf("V_0=%e\n", v_0);
	printf("STEP COUNT=%e\n", steps);
	printf("RADIAN STEP=%e\n", k);

	statsOutput << alg << ";";
	statsOutput << dp << ";";
	statsOutput << steps << ";";
	statsOutput << (double)k << ";";

	REAL alpha = mu * (pow(specificAngularMomentumSquared, -1));
	REAL u_n = 0;
	REAL beta = 3 * mu * (pow(lightningSpeedSquared, -1));

	printf("Start!\n");
	auto startTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	if (alg == 1)
	{
		if (dp)
		{
			u_n = algorithm1GR<double>(u_0, v_0, k, alpha, beta, steps, pGraph);
		}
		else
		{
			u_n = algorithm1GR<float>(u_0, v_0, k, alpha, beta, steps, pGraph);
		}
	}
	else if (alg == 2)
	{
		if (dp)
		{
			u_n = algorithm2<double>(u_0, v_0, k, alpha, steps, pGraph);
		}
		else
		{
			u_n = algorithm2<float>(u_0, v_0, k, alpha, steps, pGraph);
		}
	}
	auto finishTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	auto tiempoDeCorrida = std::chrono::duration_cast<std::chrono::nanoseconds>(finishTimer - startTimer).count();
	printf("End!\n");

	/*
	REAL ddr = 0;
	if (dp)
	{
		ddr = computeDDR<double>(u_n, semiMinorAxis, epsilon);
	}
	else
	{
		ddr = computeDDR<float>(u_n, semiMinorAxis, epsilon);
	}
	printf("DDR=%e\n", ddr);
	*/
	printf("Saving stats...\n");
	//statsOutput << ddr << ";" << (double)tiempoDeCorrida << ";" << "\n";
	statsOutput.close();

	if (print)
	{
		pGraph->EndDraw();
		pGraph->Present();

		std::wstring pngFileName = _bstr_t(runDatabaseFileName);
		pngFileName.append(L".png");
		pGraph->SavePNG(pngFileName.c_str());
	}

	printf("Done!\n");
	return 0;
}