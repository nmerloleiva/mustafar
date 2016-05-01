#include "stdafx.h"
#include "D2D1Graph.h"

using namespace std;

int main(int argc, char* argv[])
{
	double scale = 0;
	double width = 0;
	double height = 0;
	string filename;
	for (int argIndex = 1; argIndex < argc; argIndex++)
	{
		if (!strcmp(argv[argIndex], "-csv") &&
			(argIndex + 1) < argc)
		{
			filename = argv[argIndex + 1];
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

	if (scale == 0 ||
		width == 0 ||
		height == 0 ||
		filename.length() == 0)
	{
		printf("Wrong params!\n");
		printf("-csv \t Input database filename.\n");
		printf("-scale \t Drawing scale.\n");
		printf("-width \t Drawing width.\n");
		printf("-height \t Drawing height.\n");	
		return 0;
	}

	std::fstream input;
	input.open(filename, std::ios_base::in);
	if (input.fail())
	{
		printf("Error opening file!\n");
		return -1;
	}

	CD2D1Graph graph(width, height);
	graph.BeginDraw();

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

	std::wstring pngFileName = _bstr_t(filename.c_str());
	pngFileName.append(L".png");
	graph.SavePNG(pngFileName.c_str());

	printf("Done!\n");

	input.close();

	return 0;
}