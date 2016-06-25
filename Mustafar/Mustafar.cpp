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

void lagrange(double x, double X[], double y[], int Lit)
{

	double r = 0, num = 1, den = 1;
	for (int i = 0; i<Lit; i++){ //para el total de polinomios
		for (int j = 0; j<Lit; j++){ //para cada polinomio
			if (i != j){ num *= (x - X[j]); den *= (X[i] - X[j]); }
		}
		num *= y[i];
		printf("Interacion %d valor %lf\n", i, num / den);
		_getch();
		r += num / den;
		num = den = 1;
	}
	printf("\nEl resultado es: %lf", r);
}

void pointsInterpolator() //seria para ejecutar el main este proceso
{
	int m;
	double *X, *Y, x;
	void clrscr();
	printf("cuantas entradas tendra la tabla?\n\t\t");
	scanf_s("%d", &m);
	X = (double*)malloc(sizeof(double)*m);
	printf("Ingresa la tabla los valores de X:\n");
	for (int i = 0; i<m; i++) scanf_s("%lf", &X[i]);
	printf("\nIngresa la tabla los valores de Y:\n");
	Y = (double*)malloc(sizeof(double)*m);
	for (int i = 0; i<m; i++) scanf_s("%lf", &Y[i]);
	printf("Escribe el valor X para el cual se encontrara el valor de Y\n");
	scanf_s("%lf", &x);
	lagrange(x, X, Y, m);
	_getch();
}

double computeAreaError(double inicialIntervalPoint, double finalIntervalPoint, double steps){
	double error, numerator, denominator;

	numerator = pow(finalIntervalPoint - inicialIntervalPoint, 3);
	denominator = 24 * (pow(steps, 2));

	error = (numerator / denominator); // * (multiplicado por) derivada segunda evaluada en un punto del intervalo 

	return error;
}

double computePeriodError(double areaError, double steps){
	double error;
	error = (2 / steps) * areaError;  //asumo que el error del paso es 0
	return error;
}

double periodCalculation(double area, int steps){
	double period;
	period = (area * 2) / steps;
	return period;
}

double rectangularIntegral(double inicialIntervalPoint, double finalIntervalPoint, double finalIntervalPointBis, int steps) {

	double aboveArea, belowArea, abovex_i, belowx_i, aboveH, belowH, totalArea;
	inicialIntervalPoint = 0;
	finalIntervalPoint = 0, 3; // perihelio (puse un numero cualquiera)
	finalIntervalPointBis = 0, 4; // afelio (puse un numero cualquiera)
	steps = 10000000;          // Use larger value for better approximation

	aboveH = (finalIntervalPoint - inicialIntervalPoint) / steps; // Computar ancho del intervalo arriba
	belowH = (finalIntervalPointBis - finalIntervalPoint) / steps; // Computar ancho del intervalo abajo
	aboveArea = 0.0;                        // Clear running area
	belowArea = 0.0;						// Clear running area

	for (int i = 1; i <= steps - 1; i++)
	{
		abovex_i = inicialIntervalPoint + steps*aboveH;
		belowx_i = finalIntervalPoint + steps*belowH;
		aboveArea = aboveArea + (aboveH * (abovex_i));   // f(x_i) = x (abovex_i) para este ejemplo
		belowArea = belowArea + (belowH * (belowx_i));   // f(x_i) = x (belowx_i) para este ejemplo

		totalArea = aboveArea + belowArea;

		double areaErrorAbove, areaErrorBelow, totalAreaError;
		areaErrorAbove = computeAreaError(inicialIntervalPoint, finalIntervalPoint, steps);
		areaErrorBelow = computeAreaError(finalIntervalPoint, finalIntervalPointBis, steps);

		totalAreaError = areaErrorAbove + areaErrorBelow;

		// Guardar resultados
		std::fstream statsOutput;
		statsOutput.open("mustafar_solve_area_period.csv", std::ios_base::app);
		statsOutput << totalArea << ";";
		statsOutput << totalAreaError << ";";
		statsOutput << periodCalculation(totalArea, steps) << ";";
		statsOutput << computePeriodError(totalAreaError,steps) << ";";
		statsOutput << "\n";
		statsOutput.close();
	}

	return 0;
}

float f(float x) //function example
{
	return(1 / (1 + pow(x, 2)));
}

void trapezoidalIntegral()
{
	int i, n;
	float x0, xn, h, y[20], so, se, ans, x[20];
	printf("\n Enter values of x0,xn,h:\n");
	scanf_s("%f%f%f", &x0, &xn, &h);
	n = (xn - x0) / h;
	if (n % 2 == 1)
	{
		n = n + 1;
	}
	h = (xn - x0) / n;
	printf("\nrefined value of n and h are:%d  %f\n", n, h);
	printf("\n Y values \n");
	for (i = 0; i <= n; i++)
	{
		x[i] = x0 + i*h;
		y[i] = f(x[i]);
		printf("\n%f\n", y[i]);
	}
	so = 0;
	se = 0;
	for (i = 1; i<n; i++)
	{
		if (i % 2 == 1)
		{
			so = so + y[i];
		}
		else
		{
			se = se + y[i];
		}
	}
	ans = h / 3 * (y[0] + y[n] + 4 * so + 2 * se);
	printf("\nfinal integration is %f", ans);
	_getch();
}

float simpsonsOneIntegral(float(*f)(float x), float a, float b, int n) {
	float h = (b - a) / n;
	float x;
	float r;
	char m = 0;
	float s = 0.0;

	for (x = a; x <= b; x += h) {
		r = f(x);
		if (x == a || x == b) {
			s += r;
		}
		else {
			m = !m;
			s += r * (m + 1) * 2.0;
		}
	}
	return s * (h / 3.0);
}

float simpsonsTwoIntegral(float(*f)(float x), float a, float b, int n)
{
	float h = (b - a) / n;
	float sum1 = f(a + h / 2);
	float sum2 = 0;

	for (int i = 1; i < n - 1; i++)
	{
		sum1 = sum1 + f(a + h * i + h / 2);
		sum2 = sum2 + f(a + h * i);
	}
	float answer = (h / 6) * (f(a) + f(b) + 4 * sum1 + 2 * sum2);
	return answer;
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