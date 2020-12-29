#include "stdafx.h"

#include "G53IDS Bin Packing.h"

int traverseMarkovChain(int precisionLevel, int currentStateIndex, std::vector<std::vector<int>> transitionMatrix)
	{
	int newState = 0;
	int randomNumber = rand() % precisionLevel;
	int probabilitySum = 0;

	std::vector<int> transitionMatrixRow = transitionMatrix.at(currentStateIndex);

	int numberOfTransitions = transitionMatrixRow.size();

	int i = 0;
	for(i = 0; i < numberOfTransitions; i++)
		{	
		probabilitySum += transitionMatrixRow.at(i);

		if(probabilitySum >= randomNumber)
			{
			newState = i;
			break;
			}
		}

	return newState;
	}

//creates a markov model by attempting to reconstruct the transition matrix
void createMarkovModel(std::vector<int> itemStream, int precisionLevel)
	{
	//initialise transition matrix
	std::vector<std::vector<float>> transitionMatrix;

	std::vector<std::vector<int>> numberOfTransitions;

	//get unique elements from the stream of items
	std::vector<int> uniqueItemStream = itemStream;

	std::vector<int>::iterator it;
	it = std::unique(uniqueItemStream.begin(), uniqueItemStream.end());
	
	uniqueItemStream.resize(std::distance(uniqueItemStream.begin(), it));

	//sort unique item list
	std::sort(uniqueItemStream.begin(), uniqueItemStream.end());

	//vector which contains the total number of transitions from each state in the Markov Chain from each state
	std::vector<int> numberOfTransitionsFromEachState;
	
	//get difference between the highest and lowest item sizes
	int lastElementIndex = uniqueItemStream.size() - 1;
	int numberOfElements = uniqueItemStream.at(lastElementIndex) - uniqueItemStream.at(0);
	
	int firstElement = uniqueItemStream.at(0);
	int lastElement = uniqueItemStream.at(lastElementIndex);

	//fill number of transitions vector with zeroes
	int i = 0;
	for(i=0; i<(numberOfElements+1); i++)
		{
		numberOfTransitionsFromEachState.push_back(0);
		}
	
	//fill transition matrices with zeroes
	i = 0;
	for (i=0; i<(numberOfElements+1); i++)
		{
		int j = 0;
		std::vector<float> row;
		std::vector<int> row2;

		for (j=0; j<(numberOfElements+1); j++)
			{
			row.push_back(0);
			row2.push_back(0);
			}
		
		transitionMatrix.push_back(row);
		numberOfTransitions.push_back(row2);
		}

	//for each item in item stream
	i = 0;
	int itemStreamSize = itemStream.size();
	for (i = 0; i < itemStreamSize - 1; i++)
		{
		int currentItem = itemStream.at(i);
		int nextItem;

		//get index from vector
		int currentItemIndex = 0;
		int counter1 = firstElement;

		while(counter1 != currentItem)
			{
			counter1++;
			currentItemIndex++;
			}

		if (itemStream.size() > 1)
			{
			//get index from vector
			nextItem = itemStream.at(i + 1);
			int nextItemIndex = 0;
			int counter2 = firstElement;

			while(counter2 != nextItem)
				{
				counter2++;
				nextItemIndex++;
				}

			//increment number of total transitions
			int newTransitionValue = numberOfTransitionsFromEachState.at(currentItemIndex) + 1;
			numberOfTransitionsFromEachState[currentItemIndex] = newTransitionValue;

			//increment number of unique transition
			int newUniqueTransition = numberOfTransitions[currentItemIndex][nextItemIndex] + 1;
			numberOfTransitions[currentItemIndex][nextItemIndex] = newUniqueTransition;

			//update probability
			//note that all probabilities in the row have to updated
			
			std::vector<float> transitionMatrixRow = transitionMatrix.at(currentItemIndex);
			
			int j = 0;
			int transitionMatrixRowSize = transitionMatrixRow.size();
			for(j=0; j<transitionMatrixRowSize; j++)
				{
				int currentNumberOfUniqueTransitionsBetweenStates = numberOfTransitions[currentItemIndex][j];
				float tmp1 = (float)currentNumberOfUniqueTransitionsBetweenStates / (float)newTransitionValue;
				transitionMatrixRow[j] = tmp1*precisionLevel;
				}
			
			transitionMatrix[currentItemIndex] = transitionMatrixRow;
			}
		}
	}

//generate the policy matrix for the first fit bin packing algorithm
std::vector<std::vector<int>> firstFitPolicyMatrix(int binCapacity, int minElement, int maxElement)
	{
	std::vector<std::vector<int>> policyMatrix = initialisePolicyMatrix(binCapacity, minElement, maxElement);

	//return policy matrix without editing it as first fit algorithm assigns a score of 1 to all
	//active matrix indexes

	return policyMatrix;
	}

//generate the policy matrix for the best fit bin packing algorithm 
std::vector<std::vector<int>> bestFitPolicyMatrix(int binCapacity, int minElement, int maxElement)
	{
	std::vector<std::vector<int>> policyMatrix = initialisePolicyMatrix(binCapacity, minElement, maxElement);;

	int numberOfElements = (maxElement - minElement) + 1;
	int i = 0;

	for(i=0; i<numberOfElements; i++)
		{
		//assign scores to active entries
		//how scoring works is that for a given column, starting from the highest remaining capacity,
		//the score is incremented by 1 for each lower values of remaining capacity.
		//this score system follows the best fit algorithm as the item will be placed in the fullest
		//available bin, the one with the highest score for the given instance.

		int numberOfRows = policyMatrix.size();
		
		int newScore = 1;
		int startingIndex = (numberOfRows - 1) - minElement;
		int j = 0;
		
		for(j=startingIndex; j >= 0; j--)
			{
			std::vector<int> row = policyMatrix.at(j);

			int currentScore = row.at(i);

			if(currentScore == 1)
				{
				row[i] = newScore;
				policyMatrix[j] = row;
				newScore++;
				}
			}
		}

	return policyMatrix;
	} 

//generate the policy matrix for the worst fit bin packing algorithm 
std::vector<std::vector<int>> worstFitPolicyMatrix(int binCapacity, int minElement, int maxElement)
	{
	std::vector<std::vector<int>> policyMatrix = initialisePolicyMatrix(binCapacity, minElement, maxElement);;

	int numberOfElements = (maxElement - minElement) + 1;
	int i = 0;

	for (i = 0; i<numberOfElements; i++)
		{
		//assign scores to active entries
		//how scoring works is that for a given column, starting from the highest remaining capacity,
		//the score is decremented by 1 for each lower values of remaining capacity.
		//this score system follows the worst fit algorithm as the item will be placed in the emptiest
		//available bin, the one with the highest score for the given instance.

		int numberOfRows = policyMatrix.size();

		int startingIndex = (numberOfRows - 1) - minElement;
		int newScore = (startingIndex - i) + 1;
		int j = 0;

		for (j = startingIndex; j >= 0; j--)
			{
			std::vector<int> row = policyMatrix.at(j);

			int currentScore = row.at(i);

			if (currentScore == 1)
				{
				row[i] = newScore;
				policyMatrix[j] = row;
				newScore--;
				}
			}
		}

	return policyMatrix;
	}

//initialise a new policy matrix
std::vector<std::vector<int>> initialisePolicyMatrix(int binCapacity, int minElement, int maxElement)
	{
	std::vector<std::vector<int>> policyMatrix;

	//create matrix of all zeroes
	int numberOfElements = (maxElement - minElement) + 1;
	int lowestCapacity = minElement;
	int numberOfCapacities = (binCapacity - lowestCapacity) + 1;

	int i = 0; 
	int j = 0;

	for (i=lowestCapacity; i<(binCapacity+1); i++)
		{
		std::vector<int> tmp;

		for(j=0; j<numberOfElements; j++)
			{
			tmp.push_back(0);
			}

		policyMatrix.push_back(tmp);
		}
	
	//give active indexes a score of 1
	for(i=0; i<numberOfElements; i++)
		{
		int currentElement = minElement + i;

		//set starting capacity equal to 1
		int rowIndex = numberOfCapacities-1;
		
		std::vector<int> currentRow = policyMatrix.at(rowIndex);
		currentRow[i] = 1;

		policyMatrix[rowIndex] = currentRow;

		rowIndex = rowIndex - minElement;
		int currentCapacity = rowIndex + minElement;

		//while item can fit in bin, set index equal to 1
		int itemFits = 1;
		while(itemFits == 1)
			{
			if(currentCapacity < currentElement)
				{
				break;
				}
			else
				{
				currentRow = policyMatrix.at(rowIndex);
				currentRow[i] = 1;

				policyMatrix[rowIndex] = currentRow;

				rowIndex--;
				currentCapacity = rowIndex + minElement;
				}
			}
		}

	return policyMatrix;
	}

//pack items into bins using policy matrix
std::vector<std::vector<int>> packBinUsingPolicyMatrix(std::vector<std::vector<int>> policyMatrix, std::vector<std::vector<int>> pointersToBins,
	int itemToPlace, int binCapacity, int smallestElementSize)
	{		
	int i = 0;
	int numberOfBins = pointersToBins.size();
	int binFound = 0;

	//array of scores for each remaining capacity
	//note that if the bin cannot hold the required item, the score of 0 will be given
	std::vector<int> remainingCapacityScores; 
	
	for (i = 0; i < numberOfBins; i++)
		{
		std::vector<int> tempBin = pointersToBins.at(i);

		//get sum of tempBin
		int tempBinSize = tempBin.size();
		int j = 0;
		int binCapacityUsed = 0;
		for (j = 0; j < tempBinSize; j++)
			{
			binCapacityUsed += tempBin.at(j);
			}

		int remainingCapacity = binCapacity - binCapacityUsed;

		//get score of remaining capacity for the item to place from policy matrix
		int score = 0;

		if (remainingCapacity != 0 && remainingCapacity >= smallestElementSize)
			{
			int remainingCapacityIndex = remainingCapacity - smallestElementSize;
			std::vector<int> row = policyMatrix.at(remainingCapacityIndex);

			int itemIndex = itemToPlace - smallestElementSize;
			score = row.at(itemIndex);

			//add to remaining capacity scores vector
			remainingCapacityScores.push_back(score);
			}
		else 
			{
			remainingCapacityScores.push_back(0);
			}
		}

	//get score of empty bin from policy matrix
	int emptyBinIndex = policyMatrix.size()-1;
	int currentItemIndex = itemToPlace - smallestElementSize;

	std::vector<int> lastRow = policyMatrix.at(emptyBinIndex);
	int emptyBinScore = lastRow.at(currentItemIndex);

	//add empty bin score to remaining capacity scores vector
	remainingCapacityScores.push_back(emptyBinScore);

	i = 0;
	
	//find highest score in column
	int highestScore = 0;
	int highestScoreIndex = 0;
	int highestScoreAtEmptyBinIndex = 0;
	
	for (i = 0; i < numberOfBins+1; i++)
		{
		if (remainingCapacityScores.at(i) > highestScore)
			{
			if(i == numberOfBins)
				{
				highestScoreAtEmptyBinIndex = 1;
				}

			highestScore = remainingCapacityScores.at(i);
			highestScoreIndex = i;
			binFound = 1;
			}
		}

	//create a new bin if the highest score is at the empty bin index as determined by first-fit
	if(highestScoreAtEmptyBinIndex == 1)
		{
		//create new bin and add new bin to the vector of bins
		std::vector<int> newBin;
		newBin.push_back(itemToPlace);
		pointersToBins.push_back(newBin);
		}

	if (binFound == 1 && highestScoreAtEmptyBinIndex == 0)
		{
		//add element/new item to the bin
		std::vector<int> tempBin = pointersToBins.at(highestScoreIndex);
		tempBin.push_back(itemToPlace);
		pointersToBins[highestScoreIndex] = tempBin;

		//print out individual elements in the current bin (with its index number) + number of bins + items in other bins
		//printf("Added element of size %d to existing bin %d\n", itemToPlace, highestScoreIndex + 1);
		//printf("Elements in current bin (bin %d) are: ", highestScoreIndex + 1);
		//printIntVector(tempBin);

		//printf("Number of bins present are %d\n", pointersToBins.size());
		//printAllBins(pointersToBins);
		}
	else if (binFound == 0 && highestScoreAtEmptyBinIndex == 0)
		{
		//create new bin and add new bin to the vector of bins
		std::vector<int> newBin;
		newBin.push_back(itemToPlace);
		pointersToBins.push_back(newBin);

		//print out individual elements in the current bin (with its index number) + number of bins + items in other bins
		//printf("Added element of size %d to new bin %d\n", itemToPlace, numberOfBins + 1);
		//printf("Elements in current bin (bin %d) are: ", numberOfBins + 1);
		//printIntVector(newBin);

		//printf("Number of bins present are %d\n", pointersToBins.size());
		//printAllBins(pointersToBins);
		}

	return pointersToBins;
	}
	
//calculates the fitness of an individual bin, the bin capacity used
//the fitness is given as a range from 0 to 1 where 1 corresponds to the bin being full
float getBinFitness(std::vector<int> currentBin, int binCapacity)
	{
	int currentBinSize = currentBin.size();
	int binCapacityUsed = 0;
	int i = 0;

	for(i = 0; i<currentBinSize; i++)
		{
		binCapacityUsed += currentBin.at(i);
		}

	float binFitness = (float)binCapacityUsed / (float)binCapacity;
	return binFitness;
	}

//calculates the average fitness value using the fitness values from all bins
float getAverageFitnessFromAllBins(std::vector<std::vector<int>> binVector, int binCapacity)
	{
	int numberOfBins = binVector.size();
	int i = 0;
	float binFitnessSum = 0;

	for(i = 0; i<numberOfBins; i++)
		{
		binFitnessSum += getBinFitness(binVector.at(i), binCapacity);
		}

	return binFitnessSum / numberOfBins;
	}

//print out all the elements in an int vector
void printIntVector(std::vector<int> currentBin)
	{
	int i = 0;
	int binSize = currentBin.size();
	for (i = 0; i < binSize; i++)
		{
		int temp = currentBin.at(i);
		printf("%d ", temp);
		}

	printf("\n");
	}

void printAllBins(std::vector<std::vector<int>> binVector)
	{
	int numberOfBins = binVector.size();
	int i = 0;

	for (i = 0; i<numberOfBins; i++)
		{
		printf("Elements in bin %d are ", i + 1);
		printIntVector(binVector.at(i));
		}
	}

//print transition matrix
void printTransitionMatrix(std::vector<std::vector<int>> transitionMatrix)
	{
	int numberOfRows = transitionMatrix.size();
	int i = 0;
	for(i=0; i<numberOfRows; i++)
		{
		std::vector<int> currentRow = transitionMatrix.at(i);
		int numberOfCols = currentRow.size();
		int j = 0;
		for(j=0; j<numberOfCols; j++)
			{
			printf("%d ", currentRow[j]);
			}

		printf("\n");
		}

	printf("\n");
	}

//print policy matrix
void printPolicyMatrix(std::vector<std::vector<int>> policyMatrix, int binCapacity, int smallestElementSize)
	{
	int i = 0;
	int numberOfCapacities = policyMatrix.size();
	for(i=0; i<numberOfCapacities; i++)
		{
		std::vector<int> currentRow = policyMatrix.at(i);
		
		printf("Remaining Capacity %d \t", smallestElementSize + i);
		
		int j = 0;
		int numberOfItems = currentRow.size();

		for(j=0; j<numberOfItems; j++)
			{
			printf("%d ", currentRow.at(j));
			}

		printf("\n");
		}

	printf("Items in Policy Matrix: ");
	i = 0;
	std::vector<int> tmp = policyMatrix.at(0);
	int numberOfItems = tmp.size();

	for(i=0; i<numberOfItems; i++)
		{
		printf("%d ", smallestElementSize + i);
		}

	printf("\t\n\n");
	}

std::vector<int> readVariablesFromFile(const char *filePath)
	{
	std::ifstream currentFile;
	currentFile.open(filePath);

	std::string line;
	std::vector<int> variables;

	//numbers in file correspond to: 
	// 1) precision level
	// 2) starting state index
	// 3) bin capacity
	// 4) minimum item size
	// 5) maximum item size
	// 6) number of items to place
	// 7) minimum score to use in policy matrix
	// 8) maximum score to use in policy matrix
	//note that numbers 3, 4, 5 and 6 also correspond to the UBP(...) instance, this is done to make it consistent 
	//with the paper "Policy Matrix Evolution for Generation of Heuristics"

	while (getline(currentFile, line))
		{
		//if the line is not a comment in the variables file
		if (line[0] != 47 && line[1] != 47)
			{
			std::istringstream iss(line);
			int n;
			iss >> n;

			variables.push_back(n);
			}
		}

	return variables;
	}

std::vector<std::vector<int>> readTransitionMatrixFromFile(const char *filePath)
	{
	std::ifstream currentFile;
	currentFile.open(filePath);
	
	std::string line;
	std::vector<std::vector<int>> transitionMatrix;
	
	while (getline(currentFile, line))
		{
		std::vector<int> currentline;

		std::istringstream iss(line);
		int n;
		while (iss >> n)
			{
			currentline.push_back(n);
			}	

		transitionMatrix.push_back(currentline);
		}

	return transitionMatrix;
	}

//read policy matrix from file
std::vector<std::vector<int>> readPolicyMatrixFromFile(const char *filePath)
	{
	std::vector<std::vector<int>> policyMatrix;

	std::ifstream currentFile;
	currentFile.open(filePath);

	std::string line;

	while (getline(currentFile, line))
		{
		std::vector<int> currentline;

		std::istringstream iss(line);
		int n;
		while (iss >> n)
			{
			currentline.push_back(n);
			}

		currentline.erase(currentline.begin());
		policyMatrix.push_back(currentline);
		}

	policyMatrix.pop_back();
	return policyMatrix;
	}

//write policy matrix to file
void writePolicyMatrixToFile(std::vector<std::vector<int>> policyMatrix, int binCapacity, int minElement, int maxElement)
	{
	std::string filepath = "UBP("; 
	filepath += std::to_string(binCapacity);
	filepath += ",";
	filepath += std::to_string(minElement);
	filepath += ",";
	filepath += std::to_string(maxElement);
	filepath += ")-PM";
	
	//automatically add the extension ".txt" to filepath
	std::string txt = ".txt";
	filepath += txt;

	//create file of that name
	std::ofstream fileOutput;
	fileOutput.open(filepath.c_str());

	int i = 0;
	int numberOfCapacities = policyMatrix.size();
	for (i = 0; i<numberOfCapacities; i++)
		{
		std::vector<int> currentRow = policyMatrix.at(i);

		fileOutput << minElement + i << "\t"; 

		int j = 0;
		int numberOfItems = currentRow.size();

		for (j = 0; j<numberOfItems; j++)
			{
			fileOutput << currentRow.at(j) << " ";
			}

		fileOutput << std::endl;
		}

	i = 0;
	std::vector<int> tmp = policyMatrix.at(0);
	int numberOfItems = tmp.size();

	fileOutput << "\t";

	for (i = 0; i<numberOfItems; i++)
		{
		fileOutput << minElement + i << " ";
		}

	fileOutput.close();
	}

float runBinPacker(int precisionLevel, int currentState, int binCapacity, int minItemSize, int maxItemSize, int numberOfItems,
				   std::vector<std::vector<int>> policyMatrix, std::vector<std::vector<int>> transitionMatrix)
	{
	//vector containing the stream of item sizes, this is needed in order to 
	//reconstruct the Markov Chain
	std::vector<int> streamOfItems;

	//initialising a vector which contains pointers to other vectors where each
	//vector contains numbers each representing item size
	std::vector<std::vector<int>> pointersToBins;

	//a different seed is generated each time the bin packer is called
	//this guarantees that a different stream of items is generated each the time bin packer
	//is called
	srand(time(NULL));

	int i = 0;
	for (i = 0; i<numberOfItems; i++)
		{
		int currentStateIndex = currentState - minItemSize;
		currentState = traverseMarkovChain(precisionLevel, currentStateIndex, transitionMatrix) + minItemSize;
		//int numberOfItemSizes = (maxItemSize - minItemSize) + 1;
		//currentState = (rand() % numberOfItemSizes) + minItemSize;
		//printf("currentstate = %d\n", currentState);
		streamOfItems.push_back(currentState);
		int itemSize = currentState;
		//printf("-----------------------------------------------------\n");
		//printf("Item size to add %d\n", itemSize);

		int binFound = 0;

		if (i == 0)
			{
			std::vector<int> firstBin;
			firstBin.push_back(itemSize);
			pointersToBins.push_back(firstBin);

			//printf("Added element of size %d to new bin 1\n", itemSize);
			//printf("Elements in current bin (bin 1) are: ");
			//printIntVector(firstBin);

			//printf("Number of bins present are 1\n");
			//printAllBins(pointersToBins);
			}
		else
			{
			createMarkovModel(streamOfItems, precisionLevel);
			//printPolicyMatrix(policyMatrix, binCapacity, minItemSize);
			pointersToBins = packBinUsingPolicyMatrix(policyMatrix, pointersToBins, itemSize, binCapacity, minItemSize);
			}
		}

	//printf("-----------------------------------------------------\n");

	float averageBinFitness = getAverageFitnessFromAllBins(pointersToBins, binCapacity);
	//printf("Average bin fitness is %f\n", averageBinFitness);

	//printf("-----------------------------------------------------\n");

	return averageBinFitness;
	}

std::vector<std::vector<int>> normaliseTransitionMatrix(std::vector<std::vector<double>> transitionMatrix, int precisionLevel)
	{
	std::vector<std::vector<int>> newTransitionMatrix;
	int numberOfRows = transitionMatrix.size();

	int i = 0;
	for(i=0; i<numberOfRows; i++)
		{
		std::vector<double> row = transitionMatrix.at(i);
		
		//get sum of row
		int j = 0;
		double rowSum = 0;
		int colNumber = row.size();
		for(j=0; j<colNumber; j++)
			{
			rowSum += row.at(j);
			}

		//normalise row
		j = 0;
		std::vector<int> normRow;
		int normalisedRowSum = 0;
		for(j=0; j<colNumber; j++)
			{
			double currentVal = row.at(j);
			double normalisedVal = currentVal / rowSum;
			normalisedVal = normalisedVal * (double) precisionLevel;
			int normInt = (int) round(normalisedVal);
			normRow.push_back(normInt);
			normalisedRowSum += normInt;
			}

		//if sum of row is less than precision level, increase probability of remaining
		//in its current state. for example, if the current state is 5, increase the probability
		//of staying in state 5 by the remainder (precision level % sum of row).
		//this is done such that each row sums exactly to the precision level
		if(normalisedRowSum < precisionLevel)
			{
			int remainder = precisionLevel % normalisedRowSum;
			int oldProbability = normRow.at(i);
			normRow[i] = oldProbability + remainder;
			normalisedRowSum = normalisedRowSum + remainder;		
			}	
		else if(normalisedRowSum > precisionLevel)
			{
			int remainder =  normalisedRowSum % precisionLevel;
			int oldProbability = normRow.at(i);
			normRow[i] = oldProbability - remainder;
			normalisedRowSum = normalisedRowSum - remainder;
			}
	
		newTransitionMatrix.push_back(normRow);
		}

	return newTransitionMatrix;
	}

int getNumberOfActiveIndexesInPolicyMatrix(std::vector<std::vector<int>> policyMatrix)
	{
	int numberOfActiveIndexes = 0;

	int i = 0;
	int j = 0;

	int numberOfRows = policyMatrix.size();
	int numberOfCols = policyMatrix.at(0).size();

	for(i = 0; i<numberOfRows; i++)
		{
		std::vector<int> currentRow = policyMatrix.at(i);
		for(j=0; j<numberOfCols; j++)
			{
			if(currentRow.at(j) > 0)
				{
				numberOfActiveIndexes++;
				}
			}
		}

	return numberOfActiveIndexes;
	}