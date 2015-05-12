// Implementation of Genetic Algorithm to find the global minimum
// of the Schwefel function

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "GATESTER.h"	// incls the schwefel function
#include <time.h>
#include <new>


// CONSTANTS
#define POP_SIZE	100
#define RAN_NUM		((float)rand()/(RAND_MAX+1))		// a random number between 0 and 1
#define MAX_GEN		100

struct chrom_typ 
{
	float *paramvect;
	float fitness;
	
	chrom_typ(): paramvect(nullptr),fitness(0) {}
	chrom_typ(float *params, float ftns): paramvect(params), fitness(ftns) {}
};

// PROTOTYPES
float GetFitness (float* chromosome, int dim);
void printChromo (chrom_typ obj_chromosome, int dim);

// MAIN: GA implementation
int main (void)
{
	// seed the RNG
	std::srand((int)time(NULL));	

	int dim;
	std::cout << "Enter the dimension please: ";
	std::cin >> dim;

	// Initial population
	chrom_typ pop_chrom[POP_SIZE];
	
	// index
	int i_chrom (0), i_param (0);

	int gen = 0;
	// Fill 0th generation population with randomly generated chromosomes
	// paramvect:	floats in [-500, 500]
	// fitness:		0
	for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
	{
		////DEBUG
		////check if correctly initialised
		//if (i_chrom%((int)(POP_SIZE/5)) == 0)	{if (pop_chrom[i_chrom].paramvect != nullptr) {std::cout << "ERROR\n";} }
		////DEBUGEND

		pop_chrom[i_chrom].paramvect = new (std::nothrow) float[dim];
		// check for memory allocation failure
		if (pop_chrom[i_chrom].paramvect == nullptr)
		{
			std::cout << "ERROR: could not allocate memory to pop_chrom[" << i_chrom << "].paramvect\n";
			return 1;
		}

		// memory allocation successful
		// random initialisation of chromosomes
		for (i_param=0; i_param<dim; i_param++) 
		{
			pop_chrom[i_chrom].paramvect[i_param] = (1000*RAN_NUM - 500);	// scale to the domain of the Schwefel function
			
			////DEBUG
			//if (i_chrom%((int)(POP_SIZE/5)) == 0)	{std::cout << pop_chrom[i_chrom].paramvect[i_param] << ", ";}
			////DEBUGEND
		}
	}

	// Until max gen iteration is reached: assign fitness value to current population's chromosomes -> iterate (selection,crossover,mutation) to next generation
	float tot_fitness = 0.0f;
	bool isFound = false;
	while (gen<MAX_GEN) 
	{
		for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++) 
		{
			pop_chrom[i_chrom].fitness = GetFitness(pop_chrom[i_chrom].paramvect, dim);		// call fitness function and assign it to the chromosome object
			
			// Check if global minimum (0) has been found!
			if (pop_chrom[i_chrom].fitness == -1.0f)
			{
				isFound = true;
				break;
			}

			tot_fitness += pop_chrom[i_chrom].fitness;		// total sum of fitness
		}
		
		if (isFound)
		{
			// Global min found! Stop iterating
			break;
		}

		// Global min not found. Iterate!

		gen++;
	}

	// Display the best result we have
	if (isFound) 
	{
		// The global minimum (arg1, arg2, ... , arg_dim) with Schwefel fun value at #Generation
	}
	else
	{
		// maxgen ran, best approximation to globlal minimum and (...) the Schwefel fun value there

		// Best approximation
		int ind_min, maxFitness = 0;
		for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
		{
			if (pop_chrom[i_chrom].fitness > maxFitness)
			{
				maxFitness = pop_chrom[i_chrom].fitness;
				ind_min = i_chrom;
			}
		}
		std::cout << MAX_GEN << " generations\t|\t";
		printChromo(pop_chrom[ind_min], dim);
	}

	// free memory
	for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
	{
		delete [] pop_chrom[i_chrom].paramvect;
	}

	return 0;
}


// FUNCTIONS
float GetFitness (float* chromosome, int dim)
{
	// chromosome is an array of floats of length dim
	float fitness;
	fitness = schwefel(chromosome, dim);	// this fitness we want to minimise

	if (fitness != 0)
		fitness = 1/schwefel(chromosome, dim);	// this we maximise
	else
	{
		// we are at the global minimum - 0
		fitness = -1.0f;
	}

	return fitness;
}

void printChromo (chrom_typ obj_chromosome, int dim)
{
	// prints the chromosome's parameter array and Schwefel function value in format:
	// (arg1, arg2, ... , arg_dim)	|	Schwefelvalue

	std::cout << "(";
	for (int i=0; i<dim; i++)
	{
		std::cout << obj_chromosome.paramvect[i] << ", ";
	}
	std::cout << ")\t|\t" << schwefel(obj_chromosome.paramvect,dim) << "\n";

}
