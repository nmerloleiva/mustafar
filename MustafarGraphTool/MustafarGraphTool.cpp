#include "stdafx.h"
#include "D2D1Graph.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	CD2D1Graph graph(800, 800);
	graph.BeginDraw();
	
	std::fstream input;
	input.open("..\\Mustafar\\Mustafar.csv", std::ios_base::in);
	if (input.fail())
	{
		printf("Error opening file!\n");
		return -1;
	}

	float scale = 200000000;
	while (!input.eof())
	{
		string line;
		std::getline(input, line);
		vector<string> values;
		values.push_back("");
		for each (auto c in line)
		{
			if (c == ';')
			{
				values.push_back("");
			}
			else
			{
				values.back() += c;
			}
		}

		if (values.size() >= 2)
		{
			double theta_n = stod(values[0]);
			double radius_n = stod(values[1]);
			graph.DrawPointPolar(radius_n / scale, theta_n);
		}
	}

	graph.EndDraw();
	graph.Present();
	printf("Done!\n");

	return 0;
}