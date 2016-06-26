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

struct alg_data
{
	REAL u_n_min;
	REAL u_n_min_step;
	REAL u_n_max;
	REAL u_n_max_step;
	REAL u_n_4_4_N; // 2 PI
};

struct alg_params
{
	REAL np;
	REAL lambda;
	REAL steps;
	REAL alg;
	REAL dp;
	REAL print;
	REAL scale;
	REAL width;
	REAL height;
	REAL u_0;
	REAL v_0;
	REAL k;
	REAL alpha;
	REAL beta;
	REAL epsilon;
	CD2D1Graph* pGraph;
};

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
	m1 = lambda * (1.9891e+30) + experimentalDisturbance; // kg
	m2 = lambda * (3.301e+23) + experimentalDisturbance; // kg
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

#define INIT_DATA(data, u_n) \
data.u_n_min = u_0; \
data.u_n_max = u_0;	\

#define FILL_DATA_MIN_MAX(data, u_n, step) \
if (u_n < data.u_n_min) \
{ \
	data.u_n_min = u_n;	\
	data.u_n_min_step = step; \
} \
if (u_n > data.u_n_max) \
{ \
	data.u_n_max = u_n; \
	data.u_n_max_step = step; \
}

#define FILL_DATA_U_N_2_PI(data, u_n) \
data.u_n_4_4_N = u_n; \

#define PRINT(pGraph, n, k, u_n)	 \
if (pGraph)	 \
{ \
	T theta_n = n * k; \
	T radius_n = 1 / u_n; \
	pGraph->DrawPointPolar(radius_n, theta_n); \
}

template<typename T>
alg_data algorithm1(T u_0, T v_0, T k, T alpha, double steps, CD2D1Graph* pGraph = nullptr)
{
	T u_n = u_0;
	T v_n = v_0;
	T k_alpha = k * alpha;

	alg_data data;
	INIT_DATA(data, u_n);

	for (double n = 0; n < steps; n++)
	{
		FILL_DATA_MIN_MAX(data, u_n, n);
		PRINT(pGraph, n, k, u_n);

		T u_n_1 = u_n + (k * v_n);
		T v_n_1 = v_n - (k * u_n) + k_alpha;

		u_n = u_n_1;
		v_n = v_n_1;
	}

	PRINT(pGraph, steps, k, u_n);
	FILL_DATA_U_N_2_PI(data, u_n);
	return data;
}

template<typename T>
alg_data algorithm2(T u_0, T v_0, T k, T alpha, double steps, CD2D1Graph* pGraph = nullptr)
{
	T u_n = u_0;
	T v_n = v_0;

	alg_data data;
	INIT_DATA(data, u_n);

	for (double n = 0; n < steps; n++)
	{
		FILL_DATA_MIN_MAX(data, u_n, n);
		PRINT(pGraph, n, k, u_n);

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

	PRINT(pGraph, steps, k, u_n);
	FILL_DATA_U_N_2_PI(data, u_n);
	return data;
}

template<typename T>
alg_data algorithm1GR(T u_0, T v_0, T k, T alpha, T beta, double steps, CD2D1Graph* pGraph = nullptr)
{
	alg_data data;
	T u_n = u_0;
	T v_n = v_0;
	T k_alpha = k * alpha;
	T k_beta = k * beta;
	for (double n = 0; n < steps; n++)
	{
		if (pGraph)
		{
			T theta_n = n * k;
			T radius_n = 1 / u_n;
			pGraph->DrawPointPolar(radius_n, theta_n);
		}

		T u_n_1 = u_n + (k * v_n);
		T v_n_1 = v_n - (k * u_n) + k_alpha + k_beta * (pow(u_n, 2)); // u_n^2

		u_n = u_n_1;
		v_n = v_n_1;
	}

	if (pGraph)
	{
		T theta_n = steps * k;
		T radius_n = 1 / u_n;
		pGraph->DrawPointPolar(radius_n, theta_n);
	}
	data.u_n_4_4_N = u_n;
	return data;
}

template<typename T>
alg_data algorithm2GR(T u_0, T v_0, T k, T alpha, T beta, double steps, CD2D1Graph* pGraph = nullptr)
{
	alg_data data;
	T u_n = u_0;
	T v_n = v_0;
	for (double n = 0; n < steps; n++)
	{
		if (pGraph)
		{
			T theta_n = n * k;
			T radius_n = 1 / u_n;
			pGraph->DrawPointPolar(radius_n, theta_n);
		}

		T w_1 = u_n + (k * v_n) / 2;                                                                                                                                               
		T z_1 = v_n + (k * (alpha + beta * pow(u_n, 2) - u_n)) / 2;
		T w_2 = u_n + (k * z_1) / 2;
		T z_2 = v_n + (k * (alpha + beta * pow(w_1, 2) - w_1)) / 2;
		T w_3 = u_n + (k * z_2);
		T z_3 = v_n + (k * (alpha + beta * pow(w_2, 2) - w_2));

		T u_n_1 = u_n + (k * (v_n + ((2 * z_1) + (2 * z_2) + z_3))) / 6;
		T v_n_1 = v_n + (k * (
			(alpha + beta * pow(u_n, 2) - u_n) + 
			2 * (alpha + beta * pow(w_1, 2) - w_1) + 
			2 * (alpha + beta * pow(w_2, 2) - w_2) +
			(alpha + beta * pow(w_3, 2) - w_3)
			)) / 6;

		u_n = u_n_1;
		v_n = v_n_1;
	}

	if (pGraph)
	{
		T theta_n = steps * k;
		T radius_n = 1 / u_n;
		pGraph->DrawPointPolar(radius_n, theta_n);
	}
	data.u_n_4_4_N = u_n;
	return data;
}

template<typename T>
alg_data runAlgorithm(alg_params params)
{
	if (params.alg == 1)
	{
		return algorithm1<T>(params.u_0, params.v_0, params.k, params.alpha, params.steps, params.pGraph);
	}
	else if (params.alg == 2)
	{
		return algorithm2<T>(params.u_0, params.v_0, params.k, params.alpha, params.steps, params.pGraph);
	}
	else if (params.alg == 3)
	{
		return algorithm1GR<T>(params.u_0, params.v_0, params.k, params.alpha, params.beta, params.steps, params.pGraph);
	}
	else if (params.alg == 4)
	{
		return algorithm2GR<T>(params.u_0, params.v_0, params.k, params.alpha, params.beta, params.steps, params.pGraph);
	}
	return alg_data();
}

double computeAreaError(double inicialIntervalPoint, double finalIntervalPoint, double steps){
	double error, numerator, denominator;

	numerator = pow(finalIntervalPoint - inicialIntervalPoint, 3);
	denominator = 24 * (pow(steps, 2));

	error = (numerator / denominator); // * (multiplicado por) cota de derivada segunda evaluada en un punto del intervalo

	return error;
}

double computePeriodError(double areaError, double area, double specificAngularMomentumSquared, double specificAngularMomentumSquaredError){
	double error, firstTerm, secondTerm;
	
	firstTerm = (2 / sqrt(specificAngularMomentumSquared)) * areaError;
	secondTerm = ((2 * area) * (1 / specificAngularMomentumSquared)) * specificAngularMomentumSquaredError;

	error = firstTerm + secondTerm ;
	return error;
}

double periodCalculation(double area, int h){
	double period;
	period = (area * 2) / h;
	return period;
}

double rectangularIntegralAboveBelow(double inicialIntervalPoint, double finalIntervalPoint, double finalIntervalPointBis, double specificAngularMomentumSquared, int steps) {

	double aboveArea, belowArea, abovex_i, belowx_i, aboveH, belowH, totalArea;
	double specificAngularMomentumSquaredError = 0.5e-4;
	aboveH = (finalIntervalPoint - inicialIntervalPoint) / steps; // Computar ancho del intervalo arriba
	belowH = (finalIntervalPointBis - finalIntervalPoint) / steps; // Computar ancho del intervalo abajo
	aboveArea = 0.0;                        // Clear running area
	belowArea = 0.0;						// Clear running area
	double areaErrorAbove, areaErrorBelow = 0, totalAreaError;

	for (int i = 1; i <= steps - 1; i++)
	{
		abovex_i = inicialIntervalPoint + steps*aboveH;
		belowx_i = finalIntervalPoint + steps*belowH;
		aboveArea = aboveArea + (aboveH * abovex_i);   // f(x_i) = x (abovex_i) para este ejemplo
		belowArea = belowArea + (belowH * belowx_i);   // f(x_i) = x (belowx_i) para este ejemplo
		totalArea = aboveArea + belowArea;
		areaErrorAbove = computeAreaError(inicialIntervalPoint, finalIntervalPoint, steps);
		areaErrorBelow = computeAreaError(finalIntervalPoint, finalIntervalPointBis, steps);
		totalAreaError = areaErrorAbove + areaErrorBelow;
	}

	// Guardar resultados
	std::fstream statsOutput;
	statsOutput.open("mustafar_solve_area_period.csv", std::ios_base::app);
	statsOutput << totalArea << ";";
	statsOutput << totalAreaError << ";";
	statsOutput << periodCalculation(totalArea, sqrt(specificAngularMomentumSquared)) << ";";
	statsOutput << computePeriodError(totalAreaError, totalArea, specificAngularMomentumSquared, specificAngularMomentumSquaredError) << ";";
	statsOutput << "\n";
	statsOutput.close();

	return 0;
}

double rectangularIntegralGeneral(double inicialIntervalPoint, double finalIntervalPoint, double specificAngularMomentumSquared, double semiMinorAxis, double semiMayorAxis, int steps) {

	double area, x_i, h, totalArea;
	double specificAngularMomentumSquaredError = 0.5e-4;
	h = (finalIntervalPoint - inicialIntervalPoint) / steps; // Computar ancho del intervalo arriba
	area = 0.0;                        // Clear running area
	double areaError, totalAreaError;

	for (int i = 1; i <= steps - 1; i++)
	{
		x_i = inicialIntervalPoint + steps*h;
		area = area + (h * sqrt((1 - (x_i) / pow(semiMayorAxis, 2)) * pow(semiMinorAxis, 2)));   // f(x_i) = x (abovex_i) para este ejemplo
		
		totalArea = area;
		areaError = computeAreaError(inicialIntervalPoint, finalIntervalPoint, steps);
	
		totalAreaError = areaError;
	}

	// Guardar resultados
	std::fstream statsOutput;
	statsOutput.open("mustafar_solve_area_period.csv", std::ios_base::app);
	statsOutput << totalArea << ";";
	statsOutput << totalAreaError << ";";
	statsOutput << periodCalculation(totalArea, sqrt(specificAngularMomentumSquared)) << ";";
	statsOutput << computePeriodError(totalAreaError, totalArea, specificAngularMomentumSquared, specificAngularMomentumSquaredError) << ";";
	statsOutput << "\n";
	statsOutput.close();

	return 0;
}

double computeKeplerThirdLawError(double period, double semiMajorAxis, double periodError, double deltaSemiMajorAxis){
	double firstTerm, secondTerm, error;

	firstTerm = ((2 * period) / (pow(semiMajorAxis, 3))) * periodError; // (2*T / b^3) * deltaT
	secondTerm = ((3 * (pow(period, 2))) / ((pow(semiMajorAxis, 4)))) * deltaSemiMajorAxis; // (3*T^2/b^4) * deltaAxis

	error = firstTerm + secondTerm;

	return error;
}

double computeKeplerThirdLaw(double period, double semiMajorAxis){
	double quotient, numerator, denominator;

	numerator = pow(period, 2); //periodo^2
	denominator = pow(semiMajorAxis, 3);//b^3

	quotient = numerator / denominator; // cociente

	return quotient;
}

double f(double x) //function f(x) = x
{
	return(x);
}

double computeSimpsonsRuleError(double h){
	double error, numerator, denominator;

	numerator = pow(h, 5);
	denominator = 90;

	error = (numerator / denominator); // * (multiplicado por) cota de derivada cuarta evaluada en un punto del intervalo
	return error;
}

double simpsonsRule(double inicialIntervalPoint, double finalIntervalPoint, int steps)
{
	double s1 = 0, s2 = 0; //Clear running area
	double area, h;
	h = (finalIntervalPoint - inicialIntervalPoint) / steps; // Computar ancho del intervalo arriba
	if (steps % 2 == 0)
	{
		for (int i = 1; i <= steps - 1; i++)
		{
			if (i % 2 == 0)
			{
				s1 = s1 + f(inicialIntervalPoint + i*h);
			}
			else
			{
				s2 = s2 + f(inicialIntervalPoint + i*h);
			}
		}
		area = h / 3 * (f(inicialIntervalPoint) + f(finalIntervalPoint) + 4 * s2 + 2 * s1);
		printf("The value is = %f", area);
	}
	else
	{
		printf("The rule is not appliciable");
	}
	_getch();
	return area;
}

double lagrangeInterpolation(){
	{
		double x[100], y[100], quantityOfPoints, pointToInterpolate, interpolatedPoint = 0, s = 1, t = 1;
		int i, j, cutCondition = 1;
		//Respective values of the variables x and y
		quantityOfPoints = 3;
		x[0] = 0;
		x[1] = 1;
		x[2] = 2;
		y[0] = 0;
		y[1] = 1;
		y[2] = 2;
		printf("\n\n Table entered is as follows :\n\n");
		for (i = 0; i<quantityOfPoints; i++)
		{
			printf("%0.3f\t%0.3f", x[i], y[i]);
			printf("\n");
		}

		while (cutCondition == 1)
		{
			//Value of the x to find the respective value of y
			pointToInterpolate = 3.5;
			for (i = 0; i<quantityOfPoints; i++)
			{
				s = 1;
				t = 1;
				for (j = 0; j<quantityOfPoints; j++)
				{
					if (j != i)
					{
						s = s*(pointToInterpolate - x[j]);
						t = t*(x[i] - x[j]);
					}
				}
				interpolatedPoint = interpolatedPoint + ((s / t)*y[i]);
			}
			printf("\n\n The respective value of the variable y is: %f", interpolatedPoint);
			printf("\n\n Do you want to continue?\n\n Press 1 to continue and any other key to exit");
			scanf_s("%cutCondition", &cutCondition);
		}
		return interpolatedPoint;
	}
}

void solve_A_1(alg_params params, alg_data data)
{
	// Cálculo del semi eje mayor como: a =	(radio perihelio + radio afelio) / 2
	REAL semiMajorAxis = ((1 / data.u_n_min) + (1 / data.u_n_max)) / 2;
	
	// Propagación de error
	REAL deltaSemiMajorAxis = abs(1 / (2 * (pow(data.u_n_min, 2)))) * REAL_EPSILON + 
		abs(1 / (2 * (pow(data.u_n_max, 2)))) * REAL_EPSILON;

	// Cálculo del semi eje menor como: b = a * square_root(1 - e^2) 
	REAL semiMinorAxis = semiMajorAxis * sqrt((1 - pow(params.epsilon, 2)));

	// Propagación de error
	REAL deltaSemiMinorAxis = abs(sqrt((1 - pow(params.epsilon, 2)))) * deltaSemiMajorAxis +
		abs(semiMajorAxis * pow(2 * sqrt((1 - pow(params.epsilon, 2))), -1) * 2 * params.epsilon) * REAL_EPSILON;

	// Guardar resultados

	std::fstream statsOutput;
	statsOutput.open("mustafar_solve_A_1.csv", std::ios_base::app);
	statsOutput << params.alg << ";";
	statsOutput << params.steps << ";";
	statsOutput << semiMajorAxis << ";";
	statsOutput << deltaSemiMajorAxis << ";";
	statsOutput << semiMinorAxis << ";";
	statsOutput << deltaSemiMinorAxis << ";";
	statsOutput << "\n";
	statsOutput.close();
	
	if (params.pGraph)
	{
		REAL min_theta_n = data.u_n_min_step * params.k;
		REAL min_radius_n = 1 / data.u_n_min;
		params.pGraph->DrawPointPolar(min_radius_n, min_theta_n, 3, D2D1::ColorF(1, 0, 0, 1));

		REAL max_theta_n = data.u_n_max_step * params.k;
		REAL max_radius_n = 1 / data.u_n_max;
		params.pGraph->DrawPointPolar(max_radius_n, max_theta_n, 3, D2D1::ColorF(1, 0, 0, 1));
	}
}

int main(int argc, char* argv[])
{
	alg_params params;
	params.np = 92115;
	params.lambda = params.np / LAMBDA_DEN;
	params.steps = 0;
	params.alg = 0;
	params.dp = 0;
	params.print = 0;
	params.scale = 0;
	params.width = 0;
	params.height = 0;
	params.pGraph = nullptr;
	for (int argIndex = 1; argIndex < argc; argIndex++)
	{
		if (!strcmp(argv[argIndex], "-alg") &&
			(argIndex + 1) < argc)
		{
			params.alg = atol(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-steps") &&
			(argIndex + 1) < argc)
		{
			params.steps = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-dp") &&
			(argIndex + 1) < argc)
		{
			params.dp = atol(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-print") &&
			(argIndex + 1) < argc)
		{
			params.print = atol(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-width") &&
			(argIndex + 1) < argc)
		{
			params.width = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-height") &&
			(argIndex + 1) < argc)
		{
			params.height = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-scale") &&
			(argIndex + 1) < argc)
		{
			params.scale = atof(argv[argIndex + 1]);
			argIndex++;
		}
		else if (!strcmp(argv[argIndex], "-lambda") &&
			(argIndex + 1) < argc)
		{
			params.lambda = atof(argv[argIndex + 1]);
			argIndex++;
		}
	}
	
	if (params.steps == 0 ||
		params.alg < 1
		|| params.alg > 4)
	{
		printf("Wrong params!\n");
		printf("-alg \t Algorithm number {1: Euler, 2: RK4, 3: Euler GR, 4: RK4 GR}.\n");
		printf("-steps \t Algorithm step count.\n");
		printf("-dp \t Use double precision. {0, 1}\n");
		printf("-print \t Print results. {0, 1}\n");
		printf("-scale \t Print scale.\n");
		printf("-width \t Print width.\n");
		printf("-height \t Print height.\n");
		printf("-lambda \t Override lambda.\n");
		return 0;
	}

	char runDatabaseFileName[255];
	memset(runDatabaseFileName, 0, sizeof(runDatabaseFileName) / sizeof(runDatabaseFileName[0]));
	sprintf_s(runDatabaseFileName, "alg_%ld_steps_%.0f_dp_%ld_db.csv", params.alg, params.steps, params.dp);

	printf("Creating output files...\n");

	if (params.print)
	{
		params.pGraph = new CD2D1Graph(params.width, params.height, params.scale);
		params.pGraph->BeginDraw();
	}

	//std::fstream statsOutput;
	//statsOutput.open("mustafar_stats.csv", std::ios_base::app);

	REAL G = 6.673e-11; // N m^2 / kg ^2
	REAL m1 = params.lambda * (1.9891e+30); // kg
	REAL m2 = params.lambda * (3.301e+23); // kg
	REAL M = m1 + m2; // kg
	params.epsilon = pow(params.lambda, -1) * 0.2056; // -
	if (params.epsilon <= 0 || params.epsilon >= 1)
	{
		printf("Epsilon=%e out of range!\n", params.epsilon);
		return 0;
	}
	REAL semiMinorAxis = pow(params.lambda, 2) * (5.791e+10); // m
	REAL mu = G * M; // m^3 / s^2
	REAL specificAngularMomentumSquared = semiMinorAxis * mu * (1 - pow(params.epsilon, 2)); // m^4 / s^2
	params.u_0 = pow(semiMinorAxis * (1 - params.epsilon), -1);
	params.v_0 = 0;
	params.k = (2 * M_PI) / params.steps;
	REAL lightningSpeed = 3e+8;
	REAL lightningSpeedSquared = pow(lightningSpeed, 2); // c^2

	printf("ALGORITHM=%d\n", params.alg);
	printf("DOUBLE PRECISION=%d\n", params.dp);
	printf("NP=%e\n", params.np);
	printf("LAMBDA=%e\n", params.lambda);
	printf("M1=%e\n", m1);
	printf("M2=%e\n", m2);
	printf("MU=%e\n", mu);
	printf("EPSILON=%e\n", params.epsilon);
	printf("SEMI MINOR AXIS=%e\n", semiMinorAxis);
	printf("SPECIFIC ANGULAR MOMENTUM SQUARED=%e\n", specificAngularMomentumSquared);
	printf("LIGHTNING SPEED SQUARED=%e\n", lightningSpeedSquared);
	printf("U_0=%e\n", params.u_0);
	printf("V_0=%e\n", params.v_0);
	printf("STEP COUNT=%e\n", params.steps);
	printf("RADIAN STEP=%e\n", params.k);

	//statsOutput << params.alg << ";";
	//statsOutput << params.dp << ";";
	//statsOutput << params.steps << ";";
	//statsOutput << params.k << ";";

	params.alpha = mu * (pow(specificAngularMomentumSquared, -1));
	params.beta = 3 * mu * (pow(lightningSpeedSquared, -1));
	
	alg_data data;

	printf("Start!\n");

	auto startTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	if (params.dp)
	{
		data = runAlgorithm<double>(params);
	}
	else
	{
		data = runAlgorithm<float>(params);
	}

	auto finishTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	auto tiempoDeCorrida = std::chrono::duration_cast<std::chrono::nanoseconds>(finishTimer - startTimer).count();
	printf("End!\n");

	//REAL ddr = 0;
	//if (dp)
	//{
	//	ddr = computeDDR<double>(data.u_n_4_4_N, semiMinorAxis, epsilon);
	//}
	//else
	//{
	//	ddr = computeDDR<float>(data.u_n_4_4_N, semiMinorAxis, epsilon);
	//}
	//printf("DDR=%e\n", ddr);

	//printf("Saving stats...\n");
	//statsOutput << ddr << ";";
	//statsOutput << (double)tiempoDeCorrida << ";";
	//statsOutput << 1.0 / data.u_n_min << ";";
	//statsOutput << 1.0 / data.u_n_max << ";";
	//statsOutput << "\n";
	//statsOutput.close();

	solve_A_1(params, data);


	if (params.pGraph)
	{

		params.pGraph->EndDraw();
		params.pGraph->Present();

		std::wstring pngFileName = _bstr_t(runDatabaseFileName);
		pngFileName.append(L".png");
		params.pGraph->SavePNG(pngFileName.c_str());
	}

	printf("Done!\n");
	return 0;
}