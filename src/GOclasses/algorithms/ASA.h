/*
 *  ASA.h
 *  SeGMO, a Sequential Global Multiobjective Optimiser
 *
 *  Created by Dario Izzo on 5/16/08.
 *  Copyright 2008 ¿dvanced Concepts Team (European Space Agency). All rights reserved.
 *
 */

#ifndef ASA_H
#define ASA_H

#include <vector>
#include <cmath>

#include "GOproblem.h"
#include "constants.h"
#include "go_algorithm.h"
#include "individual.h"
#include "population.h"

class ASAalgorithm: public go_algorithm {
public:
	//This method initialise the SA-AN algorithm starting and final temperature setting deafult values for
	//the StartStep, the niterTemp and the niterRange. Tcoeff is evaluated accordingly.
	ASAalgorithm(int, const GOProblem &, const double &, const double &);
	Population evolve(const Population &, GOProblem &);
protected:
	size_t niterTot;
	size_t niterTemp;
	size_t niterRange;
	double T0;
	double Tcoeff;
	double StartStep;
	size_t niterOuter;
};

class mt_asa_algorithm {
	public:
		mt_asa_algorithm(int n, const GOProblem &p, const double &Ts, const double &Tf):
			m_asa(n,p,Ts,Tf) {}
		void evolve(const Population &, GOProblem &);
	private:
		struct evolver {
			void operator()(ASAalgorithm &asa, const Population &pop, GOProblem *prob) {
				asa.evolve(Population(pop),*prob);
			}
		};
		ASAalgorithm m_asa;
};

#endif
