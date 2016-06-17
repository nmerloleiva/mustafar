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
	REAL u_n_0_4_N; // 0 PI
	REAL u_n_1_4_N; // PI / 4
	REAL u_n_2_4_N; // PI / 2
	REAL u_n_3_4_N; // 3 PI / 4
	REAL u_n_4_4_N; // 2 PI
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

template<typename T>
alg_data algorithm1(T u_0, T v_0, T k, T alpha, double steps, CD2D1Graph* pGraph = nullptr)
{
	alg_data data;
	data.u_n_0_4_N = u_0;
	T u_n = u_0;
	T v_n = v_0;
	T k_alpha = k * alpha;

	for (double n = 0; n < steps; n++)
	{
		if (n == (steps / 4))
		{
			data.u_n_1_4_N = u_n;
		}
		else if (n == (steps / 2))
		{
			data.u_n_2_4_N = u_n;
		}
		else if (n == (3 * steps / 4))
		{
			data.u_n_3_4_N = u_n;
		}

		if (pGraph)
		{
			T theta_n = n * k;
			T radius_n = 1 / u_n;
			pGraph->DrawPointPolar(radius_n, theta_n);
		}

		T u_n_1 = u_n + (k * v_n);
		T v_n_1 = v_n - (k * u_n) + k_alpha;

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
alg_data algorithm2(T u_0, T v_0, T k, T alpha, double steps, CD2D1Graph* pGraph = nullptr)
{
	alg_data data;
	data.u_n_0_4_N = u_0;
	T u_n = u_0;
	T v_n = v_0;
	for (double n = 0; n < steps; n++)
	{
		if (n == (steps / 4))
		{
			data.u_n_1_4_N = u_n;
		}
		else if (n == (steps / 2))
		{
			data.u_n_2_4_N = u_n;
		}
		else if (n == (3 * steps / 4))
		{
			data.u_n_3_4_N = u_n;
		}

		if (pGraph)
		{
			T theta_n = n * k;
			T radius_n = 1 / u_n;
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
		T theta_n = steps * k;
		T radius_n = 1 / u_n;
		pGraph->DrawPointPolar(radius_n, theta_n);
	}
	data.u_n_4_4_N = u_n;
	return data;
}

template<typename T>
alg_data algorithm1GR(T u_0, T v_0, T k, T alpha, T beta, double steps, CD2D1Graph* pGraph = nullptr)
{
	alg_data data;
	data.u_n_0_4_N = u_0;
	T u_n = u_0;
	T v_n = v_0;
	T k_alpha = k * alpha;
	T k_beta = k * beta;
	for (double n = 0; n < steps; n++)
	{
		if (n == (steps / 4))
		{
			data.u_n_1_4_N = u_n;
		}
		else if (n == (steps / 2))
		{
			data.u_n_2_4_N = u_n;
		}
		else if (n == (3 * steps / 4))
		{
			data.u_n_3_4_N = u_n;
		}

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
	data.u_n_0_4_N = u_0;
	T u_n = u_0;
	T v_n = v_0;
	for (double n = 0; n < steps; n++)
	{
		if (n == (steps / 4))
		{
			data.u_n_1_4_N = u_n;
		}
		else if (n == (steps / 2))
		{
			data.u_n_2_4_N = u_n;
		}
		else if (n == (3 * steps / 4))
		{
			data.u_n_3_4_N = u_n;
		}

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

double rectangularIntegral(double(*f)(double x), double a, double b, int steps) { //function, I[a,b]
	double step = (b - a) / steps;  // width of each small rectangle
	double area = 0.0;  // signed area
	for (int i = 0; i < steps; i++) {
		area += f(a + (i + 0.5) * step) * step; // sum up each small rectangle
	}
	return area;
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

int main(int argc, char* argv[])
{
	double np = 92115;
	double lambda = np / LAMBDA_DEN;
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
		else if (!strcmp(argv[argIndex], "-lambda") &&
			(argIndex + 1) < argc)
		{
			lambda = atof(argv[argIndex + 1]);
			argIndex++;
		}
	}
	
	if (steps == 0 ||
		alg < 1
		|| alg > 4)
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

	REAL G = 6.673e-11; // N m^2 / kg ^2
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
	REAL lightningSpeed = 3e+8;
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
	printf("LIGHTNING SPEED SQUARED=%e\n", lightningSpeedSquared);
	printf("U_0=%e\n", u_0);
	printf("V_0=%e\n", v_0);
	printf("STEP COUNT=%e\n", steps);
	printf("RADIAN STEP=%e\n", k);

	statsOutput << alg << ";";
	statsOutput << dp << ";";
	statsOutput << steps << ";";
	statsOutput << (double)k << ";";

	REAL alpha = mu * (pow(specificAngularMomentumSquared, -1));
	REAL beta = 3 * mu * (pow(lightningSpeedSquared, -1));
	alg_data data;

	printf("Start!\n");
	auto startTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	if (alg == 1)
	{
		if (dp)
		{
			data = algorithm1<double>(u_0, v_0, k, alpha, steps, pGraph);
		}
		else
		{
			data = algorithm1<float>(u_0, v_0, k, alpha, steps, pGraph);
		}
	}
	else if (alg == 2)
	{
		if (dp)
		{
			data = algorithm2<double>(u_0, v_0, k, alpha, steps, pGraph);
		}
		else
		{
			data = algorithm2<float>(u_0, v_0, k, alpha, steps, pGraph);
		}
	}
	else if (alg == 3)
	{
		if (dp)
		{
			data = algorithm1GR<double>(u_0, v_0, k, alpha, beta, steps, pGraph);
		}
		else
		{
			data = algorithm1GR<float>(u_0, v_0, k, alpha, beta, steps, pGraph);
		}
	}
	else if (alg == 4)
	{
		if (dp)
		{
			data = algorithm2GR<double>(u_0, v_0, k, alpha, beta, steps, pGraph);
		}
		else
		{
			data = algorithm2GR<float>(u_0, v_0, k, alpha, beta, steps, pGraph);
		}
	}

	auto finishTimer = std::chrono::high_resolution_clock::now(); //inicio calculo del tiempo en nanosegundos
	auto tiempoDeCorrida = std::chrono::duration_cast<std::chrono::nanoseconds>(finishTimer - startTimer).count();
	printf("End!\n");

	REAL ddr = 0;
	if (dp)
	{
		ddr = computeDDR<double>(data.u_n_4_4_N, semiMinorAxis, epsilon);
	}
	else
	{
		ddr = computeDDR<float>(data.u_n_4_4_N, semiMinorAxis, epsilon);
	}
	printf("DDR=%e\n", ddr);

	printf("Saving stats...\n");
	statsOutput << ddr << ";";
	statsOutput << (double)tiempoDeCorrida << ";";
	statsOutput << 1.0 / data.u_n_0_4_N << ";";
	statsOutput << 1.0 / data.u_n_1_4_N << ";";
	statsOutput << 1.0 / data.u_n_2_4_N << ";";
	statsOutput << 1.0 / data.u_n_3_4_N << ";";
	statsOutput << "\n";
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