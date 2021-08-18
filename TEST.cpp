#include "pbPlots.hpp" 
#include "supportLib.hpp"
#include <iostream>
#include <time.h>
#include <cmath>

using namespace std;

char space[10][10] = { {'X','X','X','X','X','X','X','X','X','X'},{'X','O',' ',' ',' ',' ',' ',' ',' ','X'},
    {'X',' ',' ',' ','X','X','X','X','X','X'},{'X','X','X',' ',' ',' ',' ',' ',' ','X'},
    {'X',' ',' ',' ',' ',' ','X',' ',' ','X'},{'X',' ','X','X',' ','X','X',' ',' ','X'},
    {'X',' ',' ','X',' ',' ',' ',' ',' ','X'},{'X',' ',' ','X',' ',' ','X','X',' ','X'},
    {'X',' ',' ',' ','X','G',' ',' ',' ','X'},{'X','X','X','X','X','X','X','X','X','X'} };
const int g_col = 5, g_row = 8;
const int GENE = 20;
const int POP_SIZE = 100;
const float CROSSOVER_PROBABILITY = 0.9;
const float MUTATION_PROBABILITY = 0.1;
const int MAX_GEN = 100;
int chromosome[POP_SIZE][GENE];
double fitness[POP_SIZE];
int parent[2][GENE];
int children[2][GENE];
int survivor[POP_SIZE][GENE];
int newChromo;
int crossoverNPoint = 3;

double bestFitness = 1000.0;
double avgFitness = 0.0;
int bestChromosome[GENE];

vector<double> genBestFitness;
vector<double> genAvgFitness;
vector<double> genNum;


void initializePopulation() {
    cout << "INITIALIZE CHROMOSOME:\n";
    for (int i = 0; i < POP_SIZE; i++)
    {
        for (int j = 0; j < GENE; j++)
        {
            chromosome[i][j] = rand() % 4 + 1;
        }
    }
}

void printChromosome() {
	for (int i = 0; i < POP_SIZE; i++)
	{
		cout << "\tC" << i + 1 << ": \t";
		for (int j = 0; j < GENE; j++)
		{
			cout << chromosome[i][j] << " ";
		}
		cout << "\n";
	}
	cout << endl;
}

void evaluateChromosome() {
	cout << "\nEVALUATE CHROMOSOME:\n\n";
	int curr_col, curr_row, step;
	double dist;
	for (int i = 0; i < POP_SIZE; i++)
	{
		curr_col = 1, curr_row = 1, step = 0;
		for (int j = 0; j < GENE; j++)
		{
			if (chromosome[i][j] == 1) {
				curr_row--;
			}
			else if (chromosome[i][j] == 2) {
				curr_row++;
			}
			else if (chromosome[i][j] == 3) {
				curr_col--;
			}
			else {
				curr_col++;
			}
			if (space[curr_row][curr_col] != 'X' && space[curr_row][curr_col] != 'G' && curr_col>=0 && curr_col<10 && curr_row >= 0 && curr_row < 10) {
				step++;
			}
			else {
				if (space[curr_row][curr_col] != 'G') {
					if (chromosome[i][j] == 1) {
						curr_row++;
					}
					else if (chromosome[i][j] == 2) {
						curr_row--;
					}
					else if (chromosome[i][j] == 3) {
						curr_col++;
					}
					else {
						curr_col--;
					}
				}
				break;
			}

		}
		dist = sqrt(pow((g_col - curr_col), 2) + pow((g_row - curr_row), 2));
		fitness[i] = (step + pow(dist, 2))/100;
		cout << "\tC" << i + 1 << " \tcol = " << curr_row << " \trow = " << curr_col << " \tDist = " << dist << " \tSteps = " << step;
		cout << "  \tFV = " << fitness[i] << endl;
	}
}

void recordBestFitness() {
	for (int c = 0;c < POP_SIZE;c++) {
		if (bestFitness > fitness[c]) {
			bestFitness = fitness[c];
			for (int g = 0;g < GENE;g++) {
				bestChromosome[g] = chromosome[c][g];
			}
		}
	}
	cout << "\nBest Fitness: " << bestFitness << "\nBest Chromosome: ";
	for (int g = 0;g < GENE;g++) {
		cout << bestChromosome[g] << " ";
	}
}

void calcAvgFitness() {
	double sum = 0.0;
	for (int c = 0;c < POP_SIZE;c++) {
		sum += fitness[c];
	}
	avgFitness = sum / POP_SIZE;
	cout << "\nAverage Fitness: " << avgFitness << endl;
}

void parentSelection() {
	cout << "\n\nPARENT SELECTION:\n";
	int player1, player2, player3;
	int indexParent[2];
	for (int i = 0; i < 2; i++)
	{
		player1 = rand() % POP_SIZE;
		do {
			player2 = rand() % POP_SIZE;
			player3 = rand() % POP_SIZE;
		} while (player1 == player2 || player1 == player3 || player2 == player3);
		cout << "\n\tR" << i + 1 << "\tPlyr1: " << player1 + 1 << "\tPlyr2: " << player2 + 1 << "\tPlyr3: " << player3 + 1;
		if (fitness[player1] < fitness[player2] && fitness[player1] < fitness[player3]) {
			indexParent[i] = player1;
		}
		else if (fitness[player2] < fitness[player1] && fitness[player2] < fitness[player3]) {
			indexParent[i] = player2;
		}
		else {
			indexParent[i] = player3;
		}
		cout << "\tParent" << i + 1 << ": " << indexParent[i] + 1;
	}
	for (int i = 0; i < 2; i++) {
		cout << "\n\tP" << i + 1 << ":\t";
		for (int j = 0; j < GENE; j++) {
			parent[i][j] = chromosome[indexParent[i]][j];
			cout << parent[i][j] << " ";
		}
	}
	cout << endl;
}

void crossover() {
	cout << "\n\nCROSSOVER:\n";
	float randNo;
	int xpoint;
	for (int row = 0;row < 2;row++) {
		for (int col = 0;col < GENE;col++) {
			children[row][col] = parent[row][col];
		}
	}
	randNo = (rand() % 101) / 100.0;
	cout << "\n\tRandNo: " << randNo;
	if (randNo <= CROSSOVER_PROBABILITY) {
		for (int i = 0;i < crossoverNPoint;i++) {
			xpoint = rand() % GENE;
			cout << "\n\tXPoint: " << xpoint;
			for (int c = xpoint;c < GENE;c++) {
				children[0][c] = parent[1][c];
				children[1][c] = parent[0][c];
			}
		}
		xpoint = rand() % GENE;
		cout << "\n\tXPoint: " << xpoint;
		for (int c = xpoint;c < GENE;c++) {
			children[0][c] = parent[1][c];
			children[1][c] = parent[0][c];
		}
		for (int row = 0;row < 2;row++) {
			cout << "\n\tChildren " << row + 1 << ": ";
			for (int col = 0;col < GENE;col++) {
				cout << children[row][col] << " ";
			}
		}
	}
	else
		cout << "\n\tRandom number exceeds Crossover Probability.\n\t[CROSSOVER DOES NOT TAKE PLACE]";
	cout << endl;
}

void mutation() {
	cout << "\n\nMUTATION:\n";
	float randVal;
	int mutationBit, m_test = 0;
	for (int i = 0;i < 2;i++) {
		randVal = (rand() % 101) / 100.0;
		cout << "\n\tRandVal" << i + 1 << ": " << randVal;
		if (randVal <= MUTATION_PROBABILITY) {
			m_test = 1;
			mutationBit = rand() % (GENE - 2);
			cout << "\tMutation bit: " << mutationBit;
			children[i][mutationBit] = rand() % 4 + 1;
			children[i][mutationBit + 1] = rand() % 4 + 1;
			children[i][mutationBit + 2] = rand() % 4 + 1;
		}
		else {
			cout << "\t[Value exceeds]";
		}
		cout << "\t[Children " << i + 1 << "]";
	}
	if (m_test == 1) {
		for (int row = 0;row < 2;row++) {
			cout << "\n\tChildren " << row + 1 << ": ";
			for (int col = 0;col < GENE;col++) {
				cout << children[row][col] << " ";
			}
		}
	}
	else {
		cout << "\n\tRandom value exceeds Mutation Probability.\n\t[MUTATION DOES NOT TAKE PLACE]";
	}
	cout << endl;
}

void survivorSelection() {
	cout << "\n\nSURVIVOR SELECTON:\n";
	for (int i = 0;i < 2;i++) {
		for (int j = 0;j < GENE;j++) {
			survivor[newChromo][j] = children[i][j];
		}
		newChromo++;
	}
	for (int row = 0;row < newChromo;row++) {
		cout << "\n\tNew Chromo " << row + 1 << ":\t";
		for (int col = 0;col < GENE;col++) {
			cout << children[row][col] << " ";
		}
	}
	cout << endl;
}

void copyChromosome() {
	for (int i = 0;i < POP_SIZE;i++) {
		for (int j = 0;j < GENE;j++) {
			chromosome[i][j] = survivor[i][j];
		}
	}
}

void printResult() {
	cout << "RESULTS:\n";
	int curr_col = 1, curr_row = 1, step = 0;
	for (int i = 0;i < GENE;i++) {
		if (bestChromosome[i] == 1) {
			curr_row--;
		}
		else if (bestChromosome[i] == 2) {
			curr_row++;
		}
		else if (bestChromosome[i] == 3) {
			curr_col--;
		}
		else {
			curr_col++;
		}
		if (space[curr_row][curr_col] != 'X' && space[curr_row][curr_col] != 'G' && curr_col >= 0 && curr_col < 10 && curr_row >= 0 && curr_row < 10) {
			step++;
			space[curr_row][curr_col] = '.';
		}
		else
		{
			break;
		}
	}
	for (int i = step + 1;i < GENE;i++) {
		bestChromosome[i] = 0;
	}
	cout << "\nBest Fitness: " << bestFitness << "\nBest Chromosome: ";
	for (int g = 0;g < GENE;g++) {
		if (bestChromosome[g] != 0)
			cout << bestChromosome[g] << " ";
	}
	cout << "\n\nSpace Visualization:\n\n";
	for (int i = 0;i < 10;i++) {
		cout << "\t\t";
		for (int j = 0;j < 10;j++) {
			cout << space[i][j] << " ";
		}
		cout << endl;
	}
}

void generateGraph() {
	RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();

	ScatterPlotSeries* series = GetDefaultScatterPlotSeriesSettings();
	series->xs = &genNum;
	series->ys = &genAvgFitness;
	series->linearInterpolation = true;
	series->lineType = toVector(L"solid");
	series->lineThickness = 2;
	series->color = CreateRGBColor(1, 0, 0);

	ScatterPlotSeries* series2 = GetDefaultScatterPlotSeriesSettings();
	series2->xs = &genNum;
	series2->ys = &genBestFitness;
	series2->linearInterpolation = true;
	series2->lineType = toVector(L"solid");
	series2->lineThickness = 2;
	series2->color = CreateRGBColor(0, 0, 1);

	ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();
	settings->width = 1000;
	settings->height = 600;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
	settings->title = toVector(L"Performance");
	settings->xLabel = toVector(L"Generation");
	settings->yLabel = toVector(L"Best Fitness / Average Fitness");
	settings->scatterPlotSeries->push_back(series);
	settings->scatterPlotSeries->push_back(series2);

	DrawScatterPlotFromSettings(imageReference, settings);

	vector<double>* pngdata = ConvertToPNG(imageReference->image);
	WriteToFile(pngdata, "plot.png");
	DeleteImage(imageReference->image);
	system("plot.png");
}

int main()
{
	srand(time(0));
	initializePopulation();
	for (int gen = 0;gen < MAX_GEN;gen++) {
		int temp = 0;
		newChromo = 0;
		cout << "\n#####################################################################\nGENERATION: " << gen + 1 << "\n#####################################################################\n\n";
		printChromosome();
		evaluateChromosome();
		recordBestFitness();
		calcAvgFitness();
		genBestFitness.push_back(bestFitness);
		genAvgFitness.push_back(avgFitness);
		genNum.push_back(gen + 1);
		for (int i = 0;i < (POP_SIZE / 2);i++) {
			parentSelection();
			crossover();
			mutation();
			survivorSelection();
		}
		copyChromosome();
	}
	printResult();
	generateGraph();
}