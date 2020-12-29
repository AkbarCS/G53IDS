#pragma once

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"

#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>

int traverseMarkovChain(int precisionLevel, int currentStateIndex, std::vector<std::vector<int>> transitionMatrix);
void createMarkovModel(std::vector<int> itemStream, int precisionLevel);
std::vector<std::vector<int>> firstFitPolicyMatrix(int binCapacity, int minElement, int maxElement);
std::vector<std::vector<int>> bestFitPolicyMatrix(int binCapacity, int minElement, int maxElement);
std::vector<std::vector<int>> worstFitPolicyMatrix(int binCapacity, int minElement, int maxElement);
std::vector<std::vector<int>> initialisePolicyMatrix(int binCapacity, int minElement, int maxElement);
std::vector<std::vector<int>> packBinUsingPolicyMatrix(std::vector<std::vector<int>> policyMatrix, std::vector<std::vector<int>> pointersToBins,
													   int itemToPlace, int binCapacity, int smallestElementSize);
float getBinFitness(std::vector<int> currentBin, int binCapacity);
float getAverageFitnessFromAllBins(std::vector<std::vector<int>> binVector, int binCapacity);
void printIntVector(std::vector<int> currentBin);
void printAllBins(std::vector<std::vector<int>> binVector);
void printTransitionMatrix(std::vector<std::vector<int>> transitionMatrix);
void printPolicyMatrix(std::vector<std::vector<int>> policyMatrix, int binCapacity, int smallestElementSize);
std::vector<int> readVariablesFromFile(const char *filePath);
std::vector<std::vector<int>> readTransitionMatrixFromFile(const char *filePath);
std::vector<std::vector<int>> readPolicyMatrixFromFile(const char *filePath);
void writePolicyMatrixToFile(std::vector<std::vector<int>> policyMatrix, int binCapacity, int minElement, int maxElement);
float runBinPacker(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize,
				   int numberOfItems, std::vector<std::vector<int>> policyMatrix, std::vector<std::vector<int>> transitionMatrix);
std::vector<std::vector<int>> normaliseTransitionMatrix(std::vector<std::vector<double>> transitionMatrix, int precisionLevel);
int getNumberOfActiveIndexesInPolicyMatrix(std::vector<std::vector<int>> policyMatrix);