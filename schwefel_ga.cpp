// Implementation of Genetic Algorithm to find the global minimum
// of the Schwefel function

//////////////////////////////////////////////
//DEBUG INFO
/*
commands for debugging follow the format:

//-------------------debug TIME (9.35) DATE (130515)
	//DEBUG here
//------------------- (<, >: for begin, end)

*/
//////////////////////////////////////////////


#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "GATESTER.h"	// incls the schwefel function
#include <time.h>
#include <new>
#include <algorithm>

// CONSTANTS
#define POP_SIZE	40		// even
#define RAN_NUM		((float)rand()/(RAND_MAX))		// a random number between 0 and 1
#define MAX_GEN		30
#define PROB_X		0.4		// crossover rate
#define PROB_MUT	0.1	// mutation rate
#define EPSILON		1e-3	// precision of float in our case
#define RUN			10000		// number of GA simulations run (until global min found)
#define DIM			2		// dimension of Schwefel function
#define Nbin		100	// number of bins to store randomly generated numbers

struct chrom_typ 
{
	float* paramvect;
	float fitness;
	
	chrom_typ(): paramvect(nullptr),fitness(0.0f) {}
	chrom_typ(float *params, float ftns): paramvect(params), fitness(ftns) {}
};

// PROTOTYPES
bool floatEqual (float a, float b);
float GetFitness (float* chromosome, int dim);
void printChromo (chrom_typ obj_chromosome, int dim);
float* Roul_sel (chrom_typ* POP_CHROMO, float tot_ftns);
void crossover (float* chromo1, float* chromo2, int dim);
void mutate (float* chromosome, int dim);

// MAIN: GA implementation
int main (void)
{

	////????????????????????????????????????
	//float mylist[] = {420.9687f,420.9687f};
	//std::cout << schwefel(mylist,2);
	//return 0;
	////????????????????????????????????????

	//// ask user for the Schwefel function dimension
	//int dim;
	//std::cout << "Enter the dimension please: ";
	//std::cin >> dim;
	////

	// Seed the RNG outside run loop
	std::srand((int)time(NULL));	

	//-------------------debug 4.31 130515
	//// Check for the RNG distribution
	//int rand_bin[Nbin] = {0};	//initialise a zero-filled array as an array of bins for RNG
	//float rand_resol = (float)(1.0f/Nbin);

	//-------------------debug 9.29 130515
	//RUN GA several times
	// initialise the best chromosome from many simulations
	int best_run;
	chrom_typ best_chrom;
	best_chrom.fitness = 0.0f;

	int run = 0;
	while (run<RUN)
	{
		std::cout << "RUN# : " << run << "\n";
		//-------------------<

		//// seed the RNG
		//std::srand((int)time(NULL));	
	
		//////////////////////
		//// CRUCIAL TESTS
		//for (int i=0; i<200; i++)
		//{
		//	std::cout << RAN_NUM << ", ";
		//}
		//std::cout << (0.0f < 0.00004f) << "\n";
		//return 0;
		////
		//////////////////////

		// Initial population
		chrom_typ pop_chrom[POP_SIZE];
	
		//// Create temporary storage for next generation's population
		//chrom_typ tmp_pop[POP_SIZE];

		////-------------------debug TIME (9.35) DATE (130515)
		//std::cout << "DEBUG 9.37\n";
		////-------------------

		// index for ith chromosome, and ith gene
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

			pop_chrom[i_chrom].paramvect = new (std::nothrow) float[DIM];
			pop_chrom[i_chrom].fitness = 0.0f;

			////-------------------debug TIME (9.35) DATE (130515)
			//if (i_chrom%(POP_SIZE/5) == 0)
			//	std::cout << "DEBUG 9.40\n";
			////-------------------

			// check for memory allocation failure
			if (pop_chrom[i_chrom].paramvect == nullptr)
			{
				std::cout << "ERROR: could not allocate memory to pop_chrom[" << i_chrom << "].paramvect\n";
				return 1;
			}

			// memory allocation successful
			// random initialisation of chromosomes
			float ran_gene;		
			for (i_param=0; i_param<DIM; i_param++) 
			{
				ran_gene = (1000*RAN_NUM - 500);	// scale to the domain of the Schwefel function
				pop_chrom[i_chrom].paramvect[i_param] = ran_gene;

				////-------------------debug 4.31 130515
				//// RANDOMNESS
				//rand_bin[(int)(RAN_NUM/rand_resol)]++;		// increment the appropriate bin

				////DEBUG
				// check if the randomly assigned number does fall within [-500,500]
				if (fabs(ran_gene) > 500.0f) {std::cout << "ERROR: RNG fail [-500,500]: " << ran_gene << "\n"; return 1;}

				//if (i_chrom%((int)(POP_SIZE/5)) == 0)	{std::cout << pop_chrom[i_chrom].paramvect[i_param] << ", ";}
				////DEBUGEND
			}
		}

		////-------------------debug 4.31 130515
		//// Summary of RNG
		//for (int i=0; i<Nbin; i++)
		//{
		//	std::cout << rand_bin[i] << ", ";
		//}
		//std::cout << "\n";

		////-------------------debug TIME (9.35) DATE (130515)
		//std::cout << "DEBUG 9.41\n";
		////------------------- (<, >: for begin, end)

		// Until max gen iteration is reached: assign fitness value to current population's chromosomes -> iterate (selection,crossover,mutation) to next generation
		float tot_fitness;
		bool isFound = false;
		int ind_globalmin;

		// The iterative algorithm
		while (true) 
		{
			////-------------------debug TIME (9.58) DATE (130515)
			//std::cout << "DEBUG 9.58 generation: " << gen << "\n";
			////-------------------

			///////////////////////
			//// DEBUG MAIN
			//std::cout << "DEBUG POINT printChromo" << "\n";
			//printChromo(pop_chrom[0], DIM);
			////
			///////////////////////

			tot_fitness = 0.0f;		// reset total fitness sum to 0
		
			for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++) 
			{
				pop_chrom[i_chrom].fitness = GetFitness(pop_chrom[i_chrom].paramvect, DIM);		// call fitness function and assign it to the chromosome object
			
				////-------------------debug TIME (9.46) DATE (130515)
				//if (i_chrom%(POP_SIZE/5) == 0)
				//	std::cout << "DEBUG 9.46: " << pop_chrom[i_chrom].fitness << "\n";
				////-------------------

				// Check if global minimum (0) has been found!
				if (floatEqual(pop_chrom[i_chrom].fitness, -1.0f))
				{

					//-------------------debug TIME (9.50) DATE (130515)
					std::cout << "DEBUG 9.50: FOUND! \n";
					//-------------------

					isFound = true;
					ind_globalmin = i_chrom;
					break;
				}

				tot_fitness += pop_chrom[i_chrom].fitness;		// total sum of fitness
			}
		
			////-------------------debug TIME (9.51) DATE (130515)
			//if (gen%(MAX_GEN/5) == 0)
			//	std::cout << "DEBUG 9.51 tot_fitness: " << tot_fitness << "\n";
			////-------------------

			//-------------------debug TIME (2.37) DATE (130515)
			//if (gen%(MAX_GEN/5)==0)
			//{
			//	std::cout << "---------------------\n";
			//	std::cout << "Gen: " << gen << "\n";
			//	for (int i=0; i<POP_SIZE; i++)
			//	{
			//		if (i%(POP_SIZE/5)==0)
			//		{
			//			std::cout << "debug 2.37 fitness: " << pop_chrom[i].fitness << "\n";
			//		}
			//	}
			//}
			////-------------------


			if (isFound || (gen==MAX_GEN))
			{

				//-------------------debug TIME (9.51) DATE (130515)
				std::cout << "DEBUG 9.53 EXIT GEN: " << gen << "\n";
				//-------------------

				// Global min found OR MAX iters reached! Stop iterating
				break;
			}


			//
			// Global min not found AND MAX iters not reached... Iterate more!
			//
			// Create temporary storage for next generation's population
			chrom_typ tmp_pop[POP_SIZE];

			////-------------------debug TIME (9.56) DATE (130515)
			//std::cout << "DEBUG 9.56 init tmp_pop\n";
			////-------------------

			float *tmp_chrom1, *tmp_chrom2;	// temporary float-array chromosomes
			tmp_chrom1 = new (std::nothrow) float[DIM];
			tmp_chrom2 = new (std::nothrow) float[DIM];

			// Creation of next generation population
			// SELECTION, CROSSOVER, MUTATION
			int tmp_counter = 0;
			while (tmp_counter < POP_SIZE)
			{
				// SELECTION
				// Roulette selection of 2 parent chromosomes
				// tmp_chrom1, tmp_chrom2 can be the same pointers!
				tmp_chrom1 = Roul_sel(pop_chrom, tot_fitness);
				tmp_chrom2 = Roul_sel(pop_chrom, tot_fitness);
			
				////-------------------debug TIME (10.05) DATE (130515)
				//std::cout << "SELECTION pass\n";
				////-------------------

				////===========================<
				//// A TRICK!
				//tmp_chrom2 = tmp_chrom1;
				//std::cout << "CROSSOVER THE SAME CHROMOSOMES...\n";

				// CROSSOVER
				crossover(tmp_chrom1, tmp_chrom2, DIM);

				//std::cout << "DONE!\n";
				////
				////===========================>

				////-------------------debug TIME (10.05) DATE (130515)
				//std::cout << "XOVER pass\n";
				////-------------------

				// MUTATION
				mutate(tmp_chrom1, DIM);
				mutate(tmp_chrom2, DIM);

				////-------------------debug TIME (10.05) DATE (130515)
				//std::cout << "MUTA pass\n";
				////-------------------

				////-------------------debug TIME (10.03) DATE (130515)
				//std::cout << "DEBUG 10.03 tmp_counter: " << tmp_counter << "\n";
				////-------------------

				// storage of genetically operated daughter chromosomes into the next gen population
				tmp_pop[tmp_counter++] = chrom_typ(tmp_chrom1, 0.0f);
				tmp_pop[tmp_counter++] = chrom_typ(tmp_chrom2, 0.0f);

				//tmp_counter += 2;		// post incremented during the above storage procedure
			}

			// copy the temporary population into the main population structure
			for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
			{
				pop_chrom[i_chrom] = tmp_pop[i_chrom];
			}

			///////////////////////
			//// DEBUG MAIN
			//std::cout << "DEBUG POINT gen: " << gen << "\n";
			//std::cout << "i_chrom: " << i_chrom << "\n";
			////
			///////////////////////

			gen++;
		}

		///////////////////////
		//// DEBUG MAIN
		//std::cout << "DEBUG POINT\n";
		////
		///////////////////////

		// Display the best result we have
		if (isFound) 
		{
			// Display: Generations ran, the global minimum (arg1, arg2, ... , arg_dim), Schwefel fun value at global minimum (=0)
			std::cout << "----GLOBAL MINIMUM FOUND----\n";
			std::cout << gen << " generations\t|\t";
			printChromo(pop_chrom[ind_globalmin],DIM);

			std::cout << "\n run: " << run << "\n";
			return 0;
		}
		else
		{
			// Dislay: maxgen ran, best approximation to globlal minimum and (...) the Schwefel fun value there

			// Best approximation
			int ind_min;
			float maxFitness = 0.0f;
			for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
			{

				///////////////////////
				//// DEBUG
				//if (i_chrom%(int)(POP_SIZE/5) == 0) {std::cout << "DEBUG POINT B fitness: " << pop_chrom[i_chrom].fitness << "\n";}
				////
				///////////////////////

				////-------------------debug TIME (2.24) DATE (130515)
				//if (i_chrom%(POP_SIZE/5)==0)
				//	std::cout << "debug 2.24 fitness: " << pop_chrom[i_chrom].fitness << "\n";
				////-------------------

				if (pop_chrom[i_chrom].fitness > maxFitness)
				{
					maxFitness = pop_chrom[i_chrom].fitness;
					ind_min = i_chrom;

					if (ind_min != 0)	// the population at the exit generation is not uniform (expect for small maxgen)
					{
						//-------------------debug TIME (1.15) DATE (130515)
						std::cout << "debug 1.15 ind_min: " << ind_min << "\n";
						//-------------------
					}
				}
			}
		
			///////////////////////
			//// DEBUG
			//std::cout << "DEBUG POINT ind_min: " << ind_min << "\n";
			//std::cout << "maxFitness: " << maxFitness << "\n";
			////
			///////////////////////

			std::cout << MAX_GEN << " generations\t|\t";
			printChromo(pop_chrom[ind_min], DIM);

			// Check for fittest chromosome throughout simulation
			if (maxFitness > best_chrom.fitness) 
			{
				// update the fittest chromosome
				best_chrom = pop_chrom[ind_min];
				best_run = run;
			}
		}

		//// free memory
		//for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
		//{
		//	delete [] pop_chrom[i_chrom].paramvect;
		//}

		run++;
	//-------------------debug 9.29 130515
	//RUN GA several times
	}
	//------------------->

	// Display best chromosome from all the runs
	std::cout << "\n\n-----------\n" << "Best run: " << best_run << "\n\n";
	printChromo(best_chrom, DIM);

	return 0;
}


//////////////////////
//
// FUNCTIONS
//
//////////////////////

// floatEqual
bool floatEqual (float a, float b)
{
	return (fabs(a-b) < EPSILON);
}


// GetFitness
float GetFitness (float* chromosome, int dim)
{
	// chromosome is an array of floats of length dim
	float fitness;
	fitness = schwefel(chromosome, dim);	// this fitness we want to minimise

	if (!floatEqual(fitness, 0.0f))
	{
		fitness = ((float)1.0f/(schwefel(chromosome, dim)));	// this we maximise
	}
	else
	{
		// we are at the global minimum - 0
		fitness = -1.0f;
	}

	///////////////////////
	//// DEBUG
	//std::cout << "DEBUG POINT fitness: " << fitness << "\n";
	////
	///////////////////////
	return fitness;
}


// printChromo
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


// Roul_sel
// A Roulette selection algorithm
float* Roul_sel (chrom_typ* POP_CHROMO, float tot_ftns)
{
	int tmp_ind = 0;
	float cum_ftns = 0.0f;		// cumulative sum of fitness
	float slice_ftns = RAN_NUM*tot_ftns;	// a randomly generated float between 0 and total fitness, used as a threshold 'slice' for random sample selection
	
	///////////////////////
	//// DEBUG
	//std::cout << "DEBUG POINT tot_ftns, slice_ftns: " << tot_ftns << ", " << slice_ftns << "\n";
	////
	///////////////////////

	if ((slice_ftns > tot_ftns)||(slice_ftns < 0.0f))
	{
		// output an ERROR message
		// slice thresh is higher than the total sum of fitness
		std::cout << "ERROR: RNG failed to generate a number between 0 and 1.\n";
		return POP_CHROMO[0].paramvect;		// an error output
	}

	// check if the cum sum has exceeded the slice threshold
	do 	
	{
		cum_ftns += POP_CHROMO[tmp_ind].fitness;
		tmp_ind++;
	} while (cum_ftns < slice_ftns);
	tmp_ind--;

	if ((tmp_ind < 0)||(tmp_ind>=POP_SIZE))
	{
		std::cout << "ERROR: tmp_ind (" << tmp_ind << ") not in range\n";
		std::cout << "slice_ftns: " << slice_ftns << "\n";
	}

	return POP_CHROMO[tmp_ind].paramvect;
}


// crossover
// A simple crossover process for two chromosomes (array of floats) of length dim
void crossover (float *chromo1, float *chromo2, int dim)
{
	// Check if chromo1 and 2 are valid pointers
	if (chromo1 == nullptr)
	{
		std::cout << "ERROR: chromo1 NULL \n";
		return;
	}
	if (chromo2 == nullptr)
	{
		std::cout << "ERROR: chromo2 NULL \n";
		return;
	}

	if (chromo1 == chromo2)
	{
		// don't worry about crossing the same parents...
		//std::cout << "Same parents\n";
		return;
	}
	else
	{
		//float tmp_gene = 0.0f;
		for (int i=0; i<dim; i++)
		{
			// traverse through each gene of the chromosome and perform x-over of gene at the set rate
			if (RAN_NUM < PROB_X)
			{
				///////////////////////
				//// DEBUG
				//std::cout << "DEBUG POINT SWAP\n";
				////
				///////////////////////

				// swap the ith gene
				float tmp_gene = chromo1[i];
			
				////-------------------debug TIME (11.22) DATE (130515)
				//std::cout << "DEBUG 11.22 chrom1: " << chromo1[i] << "\n";
				////-------------------<

				chromo1[i] = chromo2[i];

				////-------------------debug TIME (11.22) DATE (130515)
				//std::cout << "DEBUG 11.22 chrom1: " << chromo1[i] << "\n";
				////

				chromo2[i] = tmp_gene;

				////-------------------debug TIME (11.22) DATE (130515)
				//std::cout << "DEBUG 11.22 chrom2: " << chromo2[i] << "\n";
				////------------------->
			}
		}
	}
}


// mutate
// A simple mutation process that randomly reassigns genes of a chromosome
void mutate (float *chromosome, int dim)
{
	for (int i=0; i<dim; i++)
	{
		// traverse through each gene and perform random reassignment at set rate
		if (RAN_NUM < PROB_MUT)
		{
			///////////////////////
			//// DEBUG
			//std::cout << "DEBUG POINT MUTATE\n";
			////
			///////////////////////

			chromosome[i] = (1000*RAN_NUM - 500);
		}
	}
}
