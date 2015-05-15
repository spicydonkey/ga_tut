//////////////////////////////////////////////////////////////////////
// Implementation of Genetic Algorithm to find the global minimum	//
// of the Schwefel function											//
//																	//
// Author:	David Shin												//
// Date:	14/05/2015												//
//////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "GATESTER.h"		// incls the schwefel function
#include <time.h>
#include <new>
#include <algorithm>

// CONSTANTS
#define POP_SIZE	100		// No. chromos in each generation (MUST be EVEN)
#define RAN_NUM		((float)rand()/(RAND_MAX))		// a random number between 0 and 1
#define MAX_GEN		1000
#define PROB_X		0.4		// crossover rate
#define PROB_MUT	0.01		// mutation rate
#define EPSILON		1e-5	// precision of float in our case
#define RUN			1		// number of GA simulations run (until global min found)
#define DIM			5		// dimension of Schwefel function
#define Nbin		100		// number of bins to store randomly generated numbers

// SIMULATION FLAGS
#define AVGSHOW		0		// 1: show average fitness of each generation	|	0: don't
#define DETAIL		0		// 1: all chroms details at every gen			|	0: don't
#define OPERSHOW	0		// 1: show genetic operation when it happens	|	0: don't

// STRUCTURES
struct chrom_typ 
{
	float* paramvect;
	float fitness;
	
	chrom_typ(): paramvect(nullptr),fitness(0.0f) {}
	chrom_typ(float *params, float ftns): paramvect(params), fitness(ftns) {}
};


// PROTOTYPES
bool	floatEqual (float a, float b);
float	GetFitness (float* chromosome, int dim);
void	printChromo (chrom_typ obj_chromosome, int dim);
float*	Roul_sel (chrom_typ* POP_CHROMO, float tot_ftns);
void	crossover (float* chromo1, float* chromo2, int dim);
void	mutate (float* chromosome, int dim);
float	avgFitness (chrom_typ* POP_CHROMO, int dim);
void	copyChromo (const float* parent, float* daughter, int length);
void	printfl (float* ptr_fl, int length);

// MAIN: GA implementation
int main (void)
{

#if (POP_SIZE%2 != 0)
	// Check if POP_SIZE is even
	std::cout << "ERROR: POP_SIZE must be an even positive integer.\n";
	return 1;
#endif


	////DEBUG: check schwefel() function in "GATESTER.h"
	//float mylist[] = {420.9687f,420.9687f};
	//std::cout << schwefel(mylist,2);
	//return 0;
	////


#ifndef DIM
	// Check if DIM is defined
	// ask user for the Schwefel function dimension
	int DIM;
	std::cout << "Enter the dimension please: ";
	std::cin >> DIM;
#endif

	// Seed the RNG outside run loop
	std::srand((int)time(NULL));	

	//-------------------debug 4.31 130515
	//// Check for the RNG distribution
	//int rand_bin[Nbin] = {0};	//initialise a zero-filled array as an array of bins for RNG
	//float rand_resol = (float)(1.0f/Nbin);

	// RUN GA multiple times
	int run = 0, best_run;
	chrom_typ best_chrom;			// the best chromo (largest fitness) from all the GA runs
	best_chrom.fitness = 0.0f;


	while (run<RUN)
	{
		std::cout << "----------------------------------------------------------------\n";
		std::cout << "RUN# : " << run << "\n";	

		//-------------------debug 4.31 130515
		//// RNG TEST
		//for (int i=0; i<200; i++)
		//{
		//	std::cout << RAN_NUM << ", ";
		//}
		//return 0;

		// Initial population
		chrom_typ pop_chrom[POP_SIZE];

		// index for ith chromosome, and ith gene
		int i_chrom (0), i_param (0);

		int gen = 0;
		// Fill 0th generation population with randomly generated chromosomes
		// paramvect:	floats in [-500, 500]
		// fitness:		0
		for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
		{

			pop_chrom[i_chrom].paramvect = new (std::nothrow) float[DIM];	// allocate sufficient memory to paramvect
			pop_chrom[i_chrom].fitness = 0.0f;

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

				// check if the randomly assigned number does fall within [-500,500]
				if (fabs(ran_gene) > 500.0f) {std::cout << "ERROR: RNG fail [-500,500]: " << ran_gene << "\n"; return 1;}
			}
		}

		////-------------------debug 4.31 130515
		//// Summary of RNG
		//for (int i=0; i<Nbin; i++)
		//{
		//	std::cout << rand_bin[i] << ", ";
		//}
		//std::cout << "\n";


		// ITERATION
		// Until max gen iteration is reached: assign fitness value to current population's chromosomes -> iterate (selection,crossover,mutation) to next generation
		float tot_fitness;
		bool isFound = false;
		int ind_globalmin;

		// fittest chromosome in this run
		int fittest_gen;

		chrom_typ fittest_this_run;

		fittest_this_run.fitness = 0.0f;
		fittest_this_run.paramvect = new (std::nothrow) float[DIM];

		// The iterative algorithm
		while (true) 
		{
			
			float tmp_ftns;
			tot_fitness = 0.0f;		// reset total fitness sum to 0
			
			for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++) 
			{
				tmp_ftns = GetFitness(pop_chrom[i_chrom].paramvect, DIM);		// call fitness function and store temporarily
				pop_chrom[i_chrom].fitness = tmp_ftns;							// assign ftns to the chromosome object

				// Check if global minimum (0) has been found!
				if (floatEqual(pop_chrom[i_chrom].fitness, -1.0f))
				{

					//-------------------debug TIME (9.50) DATE (130515)
					std::cout << "DEBUG 9.50: Global minimum found! \n";
					//-------------------

					isFound = true;
					ind_globalmin = i_chrom;
					break;
				}

				// update the fittest chromosome in this run
				if (tmp_ftns > fittest_this_run.fitness) 
				{
					for (int i=0; i<DIM; i++) 
					{
						// copy the float array
						fittest_this_run.paramvect[i] = pop_chrom[i_chrom].paramvect[i];
					}

					fittest_this_run.fitness = tmp_ftns;
					fittest_gen = gen;
				}

				tot_fitness += tmp_ftns;		// total sum of fitness
			}
		
			
#if DETAIL == 1
			// detailed output of the population in each generation
			// print all members & fitness every generation
			std::cout << "================================\n";
			std::cout << "Gen: " << gen << "\n";

			for (int i=0; i<POP_SIZE; i++)
			{
				printChromo(pop_chrom[i], DIM);
			}
#endif

			
#if AVGSHOW == 1
			// show average fitness of each generation
			std::cout << "gen <" << gen << ">  avg ftns: " << avgFitness(pop_chrom, DIM) << "\n";
#endif

			if (isFound || (gen==MAX_GEN))
			{

				//-------------------debug TIME (9.51) DATE (130515)
				std::cout << "================================================================\n";
				std::cout << "EXIT GEN: " << gen << " / " << MAX_GEN << "\n";
				//-------------------

				// Global min found OR MAX iters reached! Stop iterating
				break;
			}


			//
			// Global min not found AND MAX iters not reached... Iterate more!
			//
			// Create temporary storage for next generation's population
			chrom_typ tmp_pop[POP_SIZE];
			
			for (int i=0; i<POP_SIZE; i++) 
			{
				// Allocate sufficient memory to each object in tmp_pop[]
				tmp_pop[i].paramvect = new (std::nothrow) float[DIM];		// allocate sufficient memory to temporary population's members

				// check for memory allocation failure
				if (tmp_pop[i].paramvect == nullptr)
				{
					std::cout << "ERROR: could not allocate memory to tmp_pop[" << i << "].paramvect\n";
					return 1;
				}
			}


			const float* ptr_chromo;				// pointer to a const chromosome
			float *tmp_chrom1, *tmp_chrom2;			// temporary float-array chromosomes
			
			// tmp_chrom1 & tmp_chrom2 are pointers 
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

				// NEED TO COPY THE CHROMOSOMES (pointers!)
				ptr_chromo = Roul_sel(pop_chrom, tot_fitness);		// Roulette select a memb (float*) from pop and get ptr_chromo to point to the same addr
				copyChromo(ptr_chromo, tmp_chrom1, DIM);			// copy the genes at ptr_chromo (parent) to tmp_chrom1 (daughter) ptr

				ptr_chromo = Roul_sel(pop_chrom, tot_fitness);		// select another parent chromosome
				copyChromo(ptr_chromo, tmp_chrom2, DIM);

				
				// CROSSOVER
				crossover(tmp_chrom1, tmp_chrom2, DIM);


				// MUTATION
				mutate(tmp_chrom1, DIM);
				mutate(tmp_chrom2, DIM);

				////-------------------debug TIME (10.03) DATE (130515)
				//std::cout << "DEBUG 10.03 tmp_counter: " << tmp_counter << "\n";
				////-------------------

				// storage of genetically operated daughter chromosomes into the next gen population
				copyChromo(tmp_chrom1, tmp_pop[tmp_counter++].paramvect, DIM);		// copy the temp stored daughter chromosome into the temp_pop's correct object's member
				copyChromo(tmp_chrom2, tmp_pop[tmp_counter++].paramvect, DIM);
			}

			// COPY the temporary population into the main population structure
			for (i_chrom=0; i_chrom<POP_SIZE; i_chrom++)
			{
				copyChromo(tmp_pop[i_chrom].paramvect, pop_chrom[i_chrom].paramvect, DIM);
			}

			// Free memory allocated to tmp_pop[]
			for (int i=0; i<POP_SIZE; i++) 
			{
				//tmp_pop[i].paramvect = new (std::nothrow) float[DIM];		// allocation
				delete [] tmp_pop[i].paramvect;
			}

			gen++;
		}


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
				if (pop_chrom[i_chrom].fitness > maxFitness)
				{
					maxFitness = pop_chrom[i_chrom].fitness;
					ind_min = i_chrom;

					if (ind_min != 0)	// the population at the exit generation is not uniform (expect for small maxgen)
					{
						//-------------------debug TIME (1.15) DATE (130515)
						std::cout << "update ind_min: " << ind_min << "\n";
						//-------------------
					}
				}
			}
		

			// display the fittest member from ths run (MAY NOT BE SURVIVING)
			std::cout << "Fittest in this run from GEN <" << fittest_gen << ">:\n";
			printChromo(fittest_this_run, DIM);
			std::cout << "\n";

			// display the fittest member from the latest generation
			std::cout << "Fittest surviving: \n"; 
			printChromo(pop_chrom[ind_min], DIM);
			std::cout << "\n";

			// Check for fittest chromosome throughout simulation
			if (maxFitness > best_chrom.fitness) 
			{
				// update the fittest chromosome
				best_chrom = pop_chrom[ind_min];
				best_run = run;
			}
		}

		run++;

	//RUN GA several times
	}
	//------------------->

	// Display best chromosome from all the runs
	std::cout << "\n----------------------------------------------------------------\n\n";
	std::cout << "TOTAL NUMBER OF RUNS: " << RUN << "\n";
	std::cout << "Best run: " << best_run << "\n";
	printChromo(best_chrom, DIM);
	std::cout << "\n\n----------------------------------------------------------------\n";
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
	std::cout << ")\t|\t" << obj_chromosome.fitness << "\t|\t" << schwefel(obj_chromosome.paramvect,dim) << "\n";

}


// Roul_sel
// A Roulette selection algorithm
// BEWARE: Returns a POINTER to the paramvect of a randomly selected member in population
float* Roul_sel (chrom_typ* POP_CHROMO, float tot_ftns)
{
	int tmp_ind = 0;
	float cum_ftns = 0.0f;					// cumulative sum of fitness
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

	
	
#if OPERSHOW == 1
	// show selection behaviour
	std::cout << "SELECTED: \n";
	printChromo(POP_CHROMO[tmp_ind], DIM);
#endif

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
	
#if OPERSHOW == 1
	// show crossover behaviour
	std::cout << "CROSSOVER:\n";
	printfl(chromo1,dim);
	printfl(chromo2,dim);
	std::cout << "-----------------------------\n";
#endif

	//float tmp_gene = 0.0f;
	for (int i=0; i<dim; i++)
	{
		// traverse through each gene of the chromosome and perform x-over of gene at the set rate
		if (RAN_NUM < PROB_X)
		{
			// swap the ith gene
			float tmp_gene = chromo1[i];
			chromo1[i] = chromo2[i];
			chromo2[i] = tmp_gene;

		}
	}

#if OPERSHOW == 1
	// show crossover behaviour
	printfl(chromo1,dim);
	printfl(chromo2,dim);
#endif

}


// mutate
// A simple mutation process that randomly reassigns genes of a chromosome
void mutate (float *chromosome, int dim)
{

#if OPERSHOW == 1
	// show mutation behaviour
	std::cout << "MUTATION:\n";
	printfl(chromosome, dim);
	std::cout << "-----------------------------\n";
#endif

	for (int i=0; i<dim; i++)
	{
		// traverse through each gene and perform random reassignment at set rate
		if (RAN_NUM < PROB_MUT)
		{
			chromosome[i] = (1000*RAN_NUM - 500);
		}
	}

#if OPERSHOW == 1
	// show mutation behaviour
	printfl(chromosome, dim);
#endif

}


// avgFitness
// Calculates the average fitness of a population
float avgFitness (chrom_typ* POP_CHROMO, int dim)
{
	float sum = 0.0f;
	for (int i=0; i<POP_SIZE; i++)
	{
		sum += POP_CHROMO[i].fitness;
	}
	
	return sum/(float)POP_SIZE;
}


// copyChromo
// Copies a parent chromosome (float*)'s genes (float array values) to a daughter chrom without sharing pointer address 
void copyChromo (const float* parent, float* daughter, int length)
{
	for (int i=0; i<length; i++)
	{
		daughter[i] = parent[i];
	}
}

// printfl
// prints an array of float ( <chrom_typ>.paramvect ) from given pointer and length
void printfl (float* ptr_fl, int length)
{
	std::cout << "( ";
	for (int i=0; i<length; i++)
	{
		std::cout << ptr_fl[i] << ", ";
	}
	std::cout << ")\n";
}
