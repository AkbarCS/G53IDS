#pragma once

#include "stdafx.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <vector>

using namespace std;
// 
//  Change any of these parameters to match your needs 
//
# define POPSIZE 24
# define MAXGENS 3
# define PXOVER 1.0
# define PMUTATION 0.0175

int numberOfVariables;

//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//

struct genotype
{
	std::vector<double> gene;
	double fitness;
	double f1;
	double f2;
	std::vector<double> upper;
	std::vector<double> lower;
	double rfitness;
	double cfitness;
};

struct genotype population[POPSIZE + 1];
struct genotype newpopulation[POPSIZE + 1];

int main();

void evolveTransitionAndPolicyMatrix(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
									 std::vector<std::vector<int>> policyMatrixCompare, int numberOfRows, int numberOfTransitions, 
									 int numberOfActiveIndexesInPolicyMatrix, int minPolicyMatrixScore, int maxPolicyMatrixScore);
void evolvePolicyMatrix(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
						std::vector<std::vector<int>> transitionMatrix, int numberOfActiveIndexesInPolicyMatrix, int minPolicyMatrixScore, int maxPolicyMatrixScore);
void crossover(int &seed);
void elitist();
void evaluate(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
			  std::vector<std::vector<int>> initialPolicyMatrix, std::vector<std::vector<int>> policyMatrixCompare);
void evaluatePolicyMatrixOnly(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
							  std::vector<std::vector<int>> initialPolicyMatrix, std::vector<std::vector<int>> transitionMatrix);
int i4_uniform_ab(int a, int b, int &seed);
void initialize(int &seed, double precisionLevel, double minPolicyMatrixScore, double maxPolicyMatrixScore, int numberOfActiveIndexesInPolicyMatrix);
void keep_the_best();
void mutate(int &seed);
double r8_uniform_ab(double a, double b, int &seed);
void report(int generation);
void selector(int &seed);
void timestamp();
void Xover(int one, int two, int &seed);