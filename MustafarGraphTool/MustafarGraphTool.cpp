#include "stdafx.h"
#include "D2D1Graph.h"

using namespace std;

int main(int argc, char* argv[])
{
	CoInitialize(NULL);

	if (argc < 2)
	{
		printf("No database file specified!\n");
		return -1;
	}

	std::fstream input;
	input.open(argv[1], std::ios_base::in);
	if (input.fail())
	{
		printf("Error opening file!\n");
		return -1;
	}

	CD2D1Graph graph(800, 800);
	graph.BeginDraw();

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

	std::wstring pngFileName = _bstr_t(argv[1]);
	pngFileName.append(L".png");
	graph.SavePNG(pngFileName.c_str());

	printf("Done!\n");

	input.close();

	return 0;
}