// Mustafar.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "D2D1Graph.h"

template<typename T>
T computeDDR(T u_n, T semiMinorAxis, T epsilon)
{
	T q = semiMinorAxis * (1 - epsilon);
	T DRR = (pow(u_n, -1) / q) - 1;
	return DRR;
}

template<typename T>
T disturbeSteps(T m1,T m2,T M, T epsilon,T semiMinorAxis) //pasar parametros ByRef
{
	m1 = lambda * (1.9891e+30) + experimentalDisturbance;; // kg
	m2 = lambda * (3.301e+23) + experimentalDisturbance;; // kg
	M = m1 + m2; // kg
	semiMinorAxis = pow(lambda, 2) * (5.791e+10) + experimentalDisturbance; // m
	epsilon = pow(lambda, -1) * 0.2056 + experimentalDisturbance; // -
	if (epsilon <= 0 || epsilon >= 1)
	{
		printf("Epsilon=%e out of range!\n", epsilon);
		return 0;
	}
}

template<typename T>
T computeCp(T disturbedResult)
{
	T originalResult = 10; //poner el resultado original aca
	T relativeError = (originalResult - disturbedResult) / originalResult;
	T C_p = relativeError * (1 / experimentalDisturbance); //numero de condicion del problema
	return C_p;
}

template<typename T>
T computeStability(T resultInOnePrecision, T resultInAnotherPrecision)
{
	//utilizar la precision y los resultados de algun paso en dos precisiones diferentes y calcular
	machineUnitErrorInAnotherPrecision = 8;
	machineUnitErrorInOnePrecision = 16;
	totalMachineError = machineUnitErrorInAnotherPrecision - machineUnitErrorInOnePrecision;
	stability = (resultInOnePrecision - resultInAnotherPrecision) / (resultInOnePrecision * totalMachineError);
	return stability;
}

template<typename T>
T computeCa(T stability, T C_p)
{
	Ca = stability / C_p;
}

template<typename T>
T algorithm1(T u_0, T v_0, T k, T alpha, double steps, CD2D1Graph* pGraph = nullptr)
{
	T u_n = u_0;
	T v_n = v_0;
	T k_alpha = k * alpha;
	for (double n = 0; n < steps; n++)
	{
		if (pGraph)
		{
			REAL theta_n = n * k;
			REAL radius_n = 1 / u_n;
			pGraph->DrawPointPolar(radius_n, theta_n);
		}

		T u_n_1 = u_n + (k * v_n);
		T v_n_1 = v_n - (k * u_n) + k_alpha;

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

	printf("Start!\n");
	auto startTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	if (alg == 1)
	{
		if (dp)
		{
			u_n = algorithm1<double>(u_0, v_0, k, alpha, steps, pGraph);
		}
		else
		{
			u_n = algorithm1<float>(u_0, v_0, k, alpha, steps, pGraph);
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

	printf("Saving stats...\n");
	statsOutput << ddr << ";" << (double)tiempoDeCorrida << ";" << "\n";
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