// genetic_algorithm.cpp : Defines the entry point for the console application.
//

#define _CRT_SECURE_NO_WARNINGS

#include "stdafx.h"
#include "genetic_algorithm.h"
#include "G53IDS Bin Packing.h"

//****************************************************************************80

int main()

//****************************************************************************80
//
//  Purpose:
//
//    MAIN supervises the genetic algorithm.
//
//  Discussion:
//
//    Each generation involves selecting the best 
//    members, performing crossover & mutation and then 
//    evaluating the resulting population, until the terminating 
//    condition is satisfied   
//
//    This is a simple genetic algorithm implementation where the 
//    evaluation function takes positive values only and the 
//    fitness of an individual is the same as the value of the 
//    objective function.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    MAXGENS is the maximum number of generations.
//
//    NVARS is the number of problem variables.
//
//    PMUTATION is the probability of mutation.
//
//    POPSIZE is the population size. 
//
//    PXOVER is the probability of crossover.                          
//
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
	//file for variables to use 
	const char *variablesFilePath = "variables.txt";
	//filepath for transition Matrix to experiment on
	const char *filePath = "transitionMatrixInput.txt";
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//initialise variables
	
	//read variables (including UBP(...) instance) from file 
	std::vector<int> variables = readVariablesFromFile(variablesFilePath);

	//to adjust the precision level of the probabilities 
	int precisionLevel = variables.at(0);
	//set the starting state of the Markov Chain
	int currentState = variables.at(1);
	//bin capacity, note: bincapacity must be greater than or equal to the largest item size within the stream of items
	int binCapacity = variables.at(2);
	//minimum item size
	int minItemSize = variables.at(3);
	//maximum item size 
	int maxItemSize = variables.at(4);
	//number of items to place
	int numberOfItems = variables.at(5);
	//set policy matrix score values 
	int minPolicyMatrixScore = variables.at(6); 
	int maxPolicyMatrixScore = variables.at(7);
	//choose whether to evolve transition matrix only (set value to 1) or evolve both the policy and transition matrix (set value to 0) 
	int evolvePolicyMatrixOnly = variables.at(8);
	//select algorithm to use to generate or use the policy matrix to compare against
	int select = 0;

	std::vector<std::vector<int>> policyMatrixCompare;
	policyMatrixCompare = initialisePolicyMatrix(binCapacity, minItemSize, maxItemSize);

	if(evolvePolicyMatrixOnly == 0)
		{
		select = variables.at(9);
		if(select == 1)
			{ policyMatrixCompare = firstFitPolicyMatrix(binCapacity, minItemSize, maxItemSize); }
		else if (select == 2)
			{ policyMatrixCompare = bestFitPolicyMatrix(binCapacity, minItemSize, maxItemSize); }
		else if (select == 3)
			{ policyMatrixCompare = worstFitPolicyMatrix(binCapacity, minItemSize, maxItemSize); }
		else
			{ policyMatrixCompare = readPolicyMatrixFromFile("policyMatrixInput.txt"); }	
		}

	//printPolicyMatrix(policyMatrixCompare, binCapacity, minItemSize);

	//calculate the number of variables by calculating the number of transitions in the transition matrix and
	//the number of active indices within the policy matrix (depending on whether the policy matrix or policy matrix
	//and transition matrix pairs were being evolved)
	int numberOfItemSizes = (maxItemSize - minItemSize) + 1;
	int numberOfActiveIndexesInPolicyMatrix = getNumberOfActiveIndexesInPolicyMatrix(policyMatrixCompare) - numberOfItemSizes;
	numberOfVariables = numberOfActiveIndexesInPolicyMatrix;
	
	//read transition matrix from file - if evolvePolicyMatrixOnly == 1
	//when evolving both the transition and policy matrix, there is no need to input a transition matrix
	std::vector<std::vector<int>> transitionMatrix;
	int numberOfRows = transitionMatrix.size();
	int numberOfTransitions;
	
	if(evolvePolicyMatrixOnly == 0)
		{
		numberOfTransitions = numberOfItemSizes * numberOfItemSizes;
		numberOfVariables += numberOfTransitions;
		}
	else if(evolvePolicyMatrixOnly == 1)
		{
		transitionMatrix = readTransitionMatrixFromFile(filePath);	
		}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int numberOfActiveEntries = numberOfVariables + numberOfItemSizes;
	double optimalMutationRate = 1.0 / (double) numberOfActiveEntries;
	printf("Number of variables in chromosome %d\n", numberOfVariables);
	printf("Optimal mutation rate is 1/%d = %.5f\n", numberOfActiveEntries, optimalMutationRate);
	
	//run genetic algorithm 
	timestamp();
	cout << "\n";
	cout << "G53IDS - Policy Matrix Bin Packing\n";
	cout << "  C++ version\n";
	cout << "  Optimisation of policy matrices and transition matrix\n";

	if (numberOfVariables < 2)
		{
		cout << "\n";
		cout << "  The crossover modification will not be available,\n";
		cout << "  since it requires 2 <= NVARS.\n";
		}
		
	if (evolvePolicyMatrixOnly == 0)
		{
		evolveTransitionAndPolicyMatrix(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, policyMatrixCompare,
			           numberOfRows, numberOfTransitions, numberOfActiveIndexesInPolicyMatrix, minPolicyMatrixScore, maxPolicyMatrixScore);
		}
	else 
		{
		evolvePolicyMatrix(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, transitionMatrix, numberOfActiveIndexesInPolicyMatrix,
						   minPolicyMatrixScore, maxPolicyMatrixScore);
		}

	cout << "\n";
	cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";	
	
	//
	//  Terminate.
	//
	cout << "\n";
	cout << "SIMPLE_GA:\n";
	cout << "  Normal end of execution.\n";
	cout << "\n";
	timestamp();	

	return 0;	
} 
//****************************************************************************80

void evolveTransitionAndPolicyMatrix(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems, 
									 std::vector<std::vector<int>> policyMatrixCompare, int numberOfRows, int numberOfTransitions, 
									 int numberOfActiveIndexesInPolicyMatrix, int minPolicyMatrixScore, int maxPolicyMatrixScore)
	{
	int numberOfItemSizes = (maxItemSize - minItemSize) + 1;

	std::vector<std::vector<int>> initialPolicyMatrix = initialisePolicyMatrix(binCapacity, minItemSize, maxItemSize);

	//fill last row with value minPolicyMatrixScore+1 for normalisation
	std::vector<int> lastRow;
	int lastRowIndex = initialPolicyMatrix.size() - 1;
	lastRow.resize(numberOfItemSizes);
	std::fill(lastRow.begin(), lastRow.end(), minPolicyMatrixScore + 1);
	initialPolicyMatrix[lastRowIndex] = lastRow;
	
	//random number generation
	srand(time(NULL));
	int seed = (rand() % 10000000) + 1;

	initialize(seed, (double)precisionLevel, (double)minPolicyMatrixScore, (double)maxPolicyMatrixScore, numberOfActiveIndexesInPolicyMatrix);
	evaluate(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, initialPolicyMatrix, policyMatrixCompare);
	keep_the_best();

	int generation;
	for (generation = 0; generation < MAXGENS; generation++)
		{
		selector(seed);
		crossover(seed);
		mutate(seed);
		report(generation);
		evaluate(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, initialPolicyMatrix, policyMatrixCompare);
		elitist();
		}

	cout << "\n";
	cout << "  Best member after " << MAXGENS << " generations:\n";
	cout << "\n";

	int i;
	for (i = 0; i < numberOfVariables; i++)
		{
		cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
		}

	printf("\n");
	
	//create policy matrix
	std::vector<std::vector<int>> policyMatrix = initialPolicyMatrix;

	int numberOfColumns = (maxItemSize - minItemSize) + 1;
	int numberOfRowsInPolicyMatrix = policyMatrix.size();

	int varCounter = 0;
	int j = 0;

	for (i = 0; i<numberOfColumns; i++)
		{
		for (j = 0; j<numberOfRowsInPolicyMatrix-1; j++)
			{
			std::vector<int> currentRow = policyMatrix.at(j);

			if (currentRow.at(i) > 0)
				{
				double score = round(population[POPSIZE].gene[varCounter]);
				varCounter++;

				policyMatrix[j][i] = (int)score;
				}
			}
		}

	printf("\n");

	//create transition matrix
	std::vector<std::vector<double>> transitionMatrixUnnormalised;

	for (i = 0; i < numberOfItemSizes; i++)
		{
		std::vector<double> temp;
		int j = 0;
		for (j = 0; j < numberOfItemSizes; j++)
			{
			int probability = population[POPSIZE].gene[varCounter];
			temp.push_back(probability);
			varCounter++;
			}
			
		transitionMatrixUnnormalised.push_back(temp);
		}
		
	printf("\n");
		
	std::vector<std::vector<int>> transitionMatrixNormalised = normaliseTransitionMatrix(transitionMatrixUnnormalised, precisionLevel);
	printTransitionMatrix(transitionMatrixNormalised);

	printPolicyMatrix(policyMatrix, binCapacity, minItemSize);
	//writePolicyMatrixToFile(policyMatrix, binCapacity, minItemSize, maxItemSize);

	printf("\n");
	}

void evolvePolicyMatrix(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
						std::vector<std::vector<int>> transitionMatrix, int numberOfActiveIndexesInPolicyMatrix, int minPolicyMatrixScore, int maxPolicyMatrixScore)
	{
	int numberOfItemSizes = (maxItemSize - minItemSize) + 1;

	std::vector<std::vector<int>> initialPolicyMatrix = initialisePolicyMatrix(binCapacity, minItemSize, maxItemSize);

	//fill last row with value minPolicyMatrixScore+1 for normalisation
	std::vector<int> lastRow;
	int lastRowIndex = initialPolicyMatrix.size() - 1;
	lastRow.resize(numberOfItemSizes);
	std::fill(lastRow.begin(), lastRow.end(), minPolicyMatrixScore + 1);
	initialPolicyMatrix[lastRowIndex] = lastRow;

	//random number generation
	srand(time(NULL));
	int seed = (rand() % 10000000) + 1;
	initialize(seed, (double)precisionLevel, (double)minPolicyMatrixScore, (double)maxPolicyMatrixScore, numberOfActiveIndexesInPolicyMatrix);
	evaluatePolicyMatrixOnly(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, initialPolicyMatrix, transitionMatrix);
	keep_the_best();

	int generation;
	for (generation = 0; generation < MAXGENS; generation++)
		{
		selector(seed);
		crossover(seed);
		mutate(seed);
		report(generation);
		evaluatePolicyMatrixOnly(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, initialPolicyMatrix, transitionMatrix);
		elitist();
		}

	cout << "\n";
	cout << "  Best member after " << MAXGENS << " generations:\n";
	cout << "\n";

	int i;
	for (i = 0; i < numberOfVariables; i++)
		{
		cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
		}

	printf("\n");

	//create policy matrix
	std::vector<std::vector<int>> policyMatrix = initialPolicyMatrix;

	int numberOfColumns = (maxItemSize - minItemSize) + 1;
	int numberOfRowsInPolicyMatrix = policyMatrix.size();

	int varCounter = 0;
	int j = 0;

	for (i = 0; i<numberOfColumns; i++)
		{
		for (j = 0; j<numberOfRowsInPolicyMatrix-1; j++)
			{
			std::vector<int> currentRow = policyMatrix.at(j);

			if (currentRow.at(i) > 0)
				{
				double score = round(population[POPSIZE].gene[varCounter]);
				varCounter++;

				policyMatrix[j][i] = (int)score;
				}
			}
		}

	printf("\n");

	printPolicyMatrix(policyMatrix, binCapacity, minItemSize);
	//writePolicyMatrixToFile(policyMatrix, binCapacity, minItemSize, maxItemSize);

	printf("\n");
	}

void crossover(int &seed)

//****************************************************************************80
// 
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
	const double a = 0.0;
	const double b = 1.0;
	int mem;
	int one;
	int first = 0;
	double x;

	for (mem = 0; mem < POPSIZE; ++mem)
	{
		x = r8_uniform_ab(a, b, seed);

		if (x < PXOVER)
		{
			++first;

			if (first % 2 == 0)
			{
				Xover(one, mem, seed);
			}
			else
			{
				one = mem;
			}

		}
	}
	return;
}
//****************************************************************************80

void elitist()

//****************************************************************************80
// 
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as 
//    the last in the array. If the best member of the current 
//    generation is worse then the best member of the previous 
//    generation, the latter one would replace the worst member 
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
{
	int i;
	double best;
	int best_mem;
	double worst;
	int worst_mem;

	best = population[0].fitness;
	worst = population[0].fitness;

	for (i = 0; i < POPSIZE - 1; ++i)
	{
		if (population[i + 1].fitness < population[i].fitness)
		{

			if (best <= population[i].fitness)
			{
				best = population[i].fitness;
				best_mem = i;
			}

			if (population[i + 1].fitness <= worst)
			{
				worst = population[i + 1].fitness;
				worst_mem = i + 1;
			}

		}
		else
		{

			if (population[i].fitness <= worst)
			{
				worst = population[i].fitness;
				worst_mem = i;
			}

			if (best <= population[i + 1].fitness)
			{
				best = population[i + 1].fitness;
				best_mem = i + 1;
			}

		}

	}
	// 
	//  If the best individual from the new population is better than 
	//  the best individual from the previous population, then 
	//  copy the best from the new population; else replace the 
	//  worst individual from the current population with the 
	//  best one from the previous generation                     
	//
	if (population[POPSIZE].fitness <= best)
	{
		for (i = 0; i < numberOfVariables; i++)
		{
			population[POPSIZE].gene[i] = population[best_mem].gene[i];
			
		}
		population[POPSIZE].fitness = population[best_mem].fitness;
		population[POPSIZE].f1 = population[best_mem].f1;
		population[POPSIZE].f2 = population[best_mem].f2;
	}
	else
	{
		for (i = 0; i < numberOfVariables; i++)
		{
			population[worst_mem].gene[i] = population[POPSIZE].gene[i];
		}
		population[worst_mem].fitness = population[POPSIZE].fitness;
		population[worst_mem].f1 = population[POPSIZE].f1;
		population[worst_mem].f2 = population[POPSIZE].f2;
	}

	return;
}
//****************************************************************************80

void evaluate(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems, 
			  std::vector<std::vector<int>> initialPolicyMatrix, std::vector<std::vector<int>> policyMatrixCompare)
//****************************************************************************80
// 
//  Purpose:
//
//    EVALUATE implements the user-defined valuation function
//
//  Discussion:
//
//    Each time this is changed, the code has to be recompiled.
//    The current function is:  calls the bin packer
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
	{
	int member;
	int i;
	int j;

	for (member = 0; member < POPSIZE; member++)
		{
		//create policy matrix	
		std::vector<std::vector<int>> policyMatrixEvolved = initialPolicyMatrix;

		int numberOfColumns = (maxItemSize - minItemSize) + 1;
		int numberOfRows = initialPolicyMatrix.size();

		int varCounter = 0;

		for (i = 0; i<numberOfColumns; i++)
			{
			for (j = 0; j<numberOfRows-1; j++)
				{
				std::vector<int> currentRow = policyMatrixEvolved.at(j);
				
				//if value at the index is 0, replace value
				//replace scores column by column (item size by item size)
				if (currentRow.at(i) > 0)
					{
					double score = round(population[member].gene[varCounter]);
					varCounter++;

					policyMatrixEvolved[j][i] = (int)score;
					}
				}
			}
		
		std::vector<std::vector<double>> transitionMatrix;
		
		//create transition matrix
		int numberOfItemSizes = (maxItemSize - minItemSize) + 1;
		for (i = 0; i < numberOfItemSizes; i++)
			{
			std::vector<double> temp;
			for (j = 0; j < numberOfItemSizes; j++)
				{
				double probability = population[member].gene[varCounter];
				temp.push_back(probability);
				varCounter++;
				}

			transitionMatrix.push_back(temp);
			}

		//normalise transition matrix
		std::vector<std::vector<int>> normTransitionMatrix;
		normTransitionMatrix = normaliseTransitionMatrix(transitionMatrix, precisionLevel);

		//calculate fitness
		double f1 = runBinPacker(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, policyMatrixEvolved, normTransitionMatrix);
		double f2 = runBinPacker(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, policyMatrixCompare, normTransitionMatrix);

		double fitness = f1 - f2;

		population[member].fitness = fitness;
		population[member].f1 = f1;
		population[member].f2 = f2;
		}
	}

void evaluatePolicyMatrixOnly(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
							  std::vector<std::vector<int>> initialPolicyMatrix, std::vector<std::vector<int>> transitionMatrix)
	{
	int member;
	int i;
	int j;

	for (member = 0; member < POPSIZE; member++)
		{
		//create policy matrix	
		std::vector<std::vector<int>> policyMatrixEvolved = initialPolicyMatrix;

		int numberOfColumns = (maxItemSize - minItemSize) + 1;
		int numberOfRows = initialPolicyMatrix.size();

		int varCounter = 0;

		for (i = 0; i<numberOfColumns; i++)
			{
			for (j = 0; j<numberOfRows-1; j++)
				{
				std::vector<int> currentRow = policyMatrixEvolved.at(j);

				//if value at the index is 0, replace value
				//replace scores column by column (item size by item size)
				if (currentRow.at(i) > 0)
					{
					double score = round(population[member].gene[varCounter]);
					varCounter++;

					policyMatrixEvolved[j][i] = (int)score;
					}
				}
			}

		//calculate fitness
		double fitness = runBinPacker(precisionLevel, currentState, binCapacity, minItemSize, maxItemSize, numberOfItems, policyMatrixEvolved, transitionMatrix);

		population[member].fitness = fitness;
		}
	}


//****************************************************************************80

int i4_uniform_ab(int a, int b, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
	int c;
	const int i4_huge = 2147483647;
	int k;
	float r;
	int value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "I4_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}
	//
	//  Guarantee A <= B.
	//
	if (b < a)
	{
		c = a;
		a = b;
		b = c;
	}

	k = seed / 127773;

	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	r = (float)(seed) * 4.656612875E-10;
	//
	//  Scale R to lie between A-0.5 and B+0.5.
	//
	r = (1.0 - r) * ((float)a - 0.5)
		+ r   * ((float)b + 0.5);
	//
	//  Use rounding to convert R to an integer between A and B.
	//
	value = round(r);
	//
	//  Guarantee A <= VALUE <= B.
	//
	if (value < a)
	{
		value = a;
	}
	if (b < value)
	{
		value = b;
	}

	return value;
}
//****************************************************************************80

void initialize(int &seed, double precisionLevel, double minPolicyMatrixScore, double maxPolicyMatrixScore, int numberOfActiveIndexesInPolicyMatrix)

//****************************************************************************80
// 
//  Purpose:
//
//    INITIALIZE initializes the genes within the variables bounds. 
//
//  Discussion:
//
//    It also initializes (to zero) all fitness values for each
//    member of the population. It reads upper and lower bounds 
//    of each variable from the input file `gadata.txt'. It 
//    randomly generates values between these bounds for each 
//    gene of each genotype in the population. The format of 
//    the input file `gadata.txt' is 
//
//      var1_lower_bound var1_upper bound
//      var2_lower_bound var2_upper bound ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
	int i;
	int j;

	for(i=0; i<POPSIZE+1; i++)
		{
		std::vector<double> genetmp;
		genetmp.resize(numberOfVariables);
		population[i].gene = genetmp;

		std::vector<double> uppertmp;
		uppertmp.resize(numberOfVariables);
		population[i].upper = uppertmp;

		std::vector<double> lowertmp;
		lowertmp.resize(numberOfVariables);
		population[i].lower = lowertmp;
		}

	// 
	//  Initialize variables within the bounds 
	//

	//initialise policy matrix
	for (i = 0; i < numberOfActiveIndexesInPolicyMatrix; i++)
		{
		for (j = 0; j < POPSIZE; j++)
			{
			population[j].fitness = 0;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			population[j].f1 = 0;
			population[j].f2 = 0;
			population[j].lower[i] = minPolicyMatrixScore;
			population[j].upper[i] = maxPolicyMatrixScore;
			population[j].gene[i] = r8_uniform_ab(minPolicyMatrixScore, maxPolicyMatrixScore, seed);
			}
		}

	//initialise transition matrix
	for (i = numberOfActiveIndexesInPolicyMatrix; i < numberOfVariables; i++)
		{
		for (j = 0; j < POPSIZE; j++)
			{
			population[j].fitness = 0;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			population[j].lower[i] = 0.0;
			population[j].upper[i] = precisionLevel;
			population[j].gene[i] = r8_uniform_ab(0.0, precisionLevel, seed);
			}
		}

	return;
}
//****************************************************************************80

void keep_the_best()

//****************************************************************************80
// 
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population. 
//
//  Discussion:
//
//    Note that the last entry in the array Population holds a 
//    copy of the best individual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int CUR_BEST, the index of the best individual.
//
{
	int cur_best;
	int mem;
	int i;

	cur_best = 0;

	for (mem = 0; mem < POPSIZE; mem++)
	{
		if (population[POPSIZE].fitness < population[mem].fitness)
		{
			cur_best = mem;
			population[POPSIZE].fitness = population[mem].fitness;
		}
	}
	// 
	//  Once the best member in the population is found, copy the genes.
	//
	for (i = 0; i < numberOfVariables; i++)
	{
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
	}

	population[POPSIZE].f1 = population[cur_best].f1;
	population[POPSIZE].f2 = population[cur_best].f2;

	return;
}
//****************************************************************************80

void mutate(int &seed)

//****************************************************************************80
// 
//  Purpose:
//
//    MUTATE performs a random uniform mutation. 
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value 
//    between the lower and upper bounds of this variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	double lbound;
	double ubound;
	double x;

	for (i = 0; i < POPSIZE; i++)
	{
		for (j = 0; j < numberOfVariables; j++)
		{
			x = r8_uniform_ab(a, b, seed);
			if (x < PMUTATION)
			{
				lbound = population[i].lower[j];
				ubound = population[i].upper[j];
				population[i].gene[j] = r8_uniform_ab(lbound, ubound, seed);
			}
		}
	}

	return;
}
//****************************************************************************80

double r8_uniform_ab(double a, double b, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_AB returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_AB, a number strictly between A and B.
//
{
	int i4_huge = 2147483647;
	int k;
	double value;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "R8_UNIFORM_AB - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	k = seed / 127773;

	seed = 16807 * (seed - k * 127773) - k * 2836;

	if (seed < 0)
	{
		seed = seed + i4_huge;
	}

	value = (double)(seed) * 4.656612875E-10;

	value = a + (b - a) * value;

	return value;
}
//****************************************************************************80

void report(int generation)

//****************************************************************************80
// 
//  Purpose:
//
//    REPORT reports progress of the simulation. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double avg, the average population fitness.
//
//    Local, best_val, the best population fitness.
//
//    Local, double square_sum, square of sum for std calc.
//
//    Local, double stddev, standard deviation of population fitness.
//
//    Local, double sum, the total population fitness.
//
//    Local, double sum_square, sum of squares for std calc.
//
{
	double avg;
	double best_val;
	int i;
	double square_sum;
	double stddev;
	double sum;
	double f1;
	double f2;
	double sum_square;

	if (generation == 0)
	{
		cout << "\n";
		cout << "  Generation       Best            Average       Standard       f1       f2\n";
		cout << "  number           value           fitness       deviation      value    value\n";
		cout << "\n";
	}

	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < POPSIZE; i++)
	{
		sum = sum + population[i].fitness;
		sum_square = sum_square + population[i].fitness * population[i].fitness;
	}

	avg = sum / (double)POPSIZE;
	square_sum = avg * avg * POPSIZE;
	stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
	best_val = population[POPSIZE].fitness;
	f1 = population[POPSIZE].f1;
	f2 = population[POPSIZE].f2;

	cout << "  " << setw(8) << generation
		<< "  " << setw(14) << best_val
		<< "  " << setw(14) << avg
		<< "  " << setw(14) << stddev
		<< "  " << setw(14) << f1
		<< "  " << setw(14) << f2 << "\n";

	return;
}
//****************************************************************************80

void selector(int &seed)

//****************************************************************************80
// 
//  Purpose:
//
//    SELECTOR is the selection function.
//
//  Discussion:
//
//    Standard proportional selection for maximization problems incorporating 
//    the elitist model.  This makes sure that the best member always survives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
	const double a = 0.0;
	const double b = 1.0;
	int i;
	int j;
	int mem;
	double p;
	double sum;
	//
	//  Find the total fitness of the population.
	//
	sum = 0.0;
	for (mem = 0; mem < POPSIZE; mem++)
	{
		sum = sum + population[mem].fitness;
	}
	//
	//  Calculate the relative fitness of each member.
	//
	for (mem = 0; mem < POPSIZE; mem++)
	{
		population[mem].rfitness = population[mem].fitness / sum;
	}
	// 
	//  Calculate the cumulative fitness.
	//
	population[0].cfitness = population[0].rfitness;
	for (mem = 1; mem < POPSIZE; mem++)
	{
		population[mem].cfitness = population[mem - 1].cfitness +
			population[mem].rfitness;
	}
	// 
	//  Select survivors using cumulative fitness. 
	//
	for (i = 0; i < POPSIZE; i++)
	{
		p = r8_uniform_ab(a, b, seed);
		if (p < population[0].cfitness)
		{
			newpopulation[i] = population[0];
		}
		else
		{
			for (j = 0; j < POPSIZE; j++)
			{
				if (population[j].cfitness <= p && p < population[j + 1].cfitness)
				{
					newpopulation[i] = population[j + 1];
				}
			}
		}
	}
	// 
	//  Overwrite the old population with the new one.
	//
	for (i = 0; i < POPSIZE; i++)
	{
		population[i] = newpopulation[i];
	}

	return;
}
//****************************************************************************80

void timestamp()

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2003
//
//  Author:
//
//    John Burkardt
//
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time(NULL);
	tm = localtime(&now);

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

	cout << time_buffer << "\n";

	return;
# undef TIME_SIZE
}
//****************************************************************************80

void Xover(int one, int two, int &seed)

//****************************************************************************80
// 
//  Purpose:
//
//    XOVER performs crossover of the two selected parents. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
{
	int i;
	int point;
	double t;
	// 
	//  Select the crossover point.
	//
	point = i4_uniform_ab(0, numberOfVariables - 1, seed);
	//
	//  Swap genes in positions 0 through POINT-1.
	//
	for (i = 0; i < point; i++)
	{
		t = population[one].gene[i];
		population[one].gene[i] = population[two].gene[i];
		population[two].gene[i] = t;
	}

	return;
}