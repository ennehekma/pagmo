/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef PAGMO_ALGORITHM_BASE_NLOPT_H
#define PAGMO_ALGORITHM_BASE_NLOPT_H

#include <cstddef>
#include <nlopt.h>
#include <string>

#include "../config.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// Base class for wrapping NLopt's algorithms.
class __PAGMO_VISIBLE base_nlopt: public base
{
	protected:
		base_nlopt(nlopt_algorithm, bool, int, const double &);
		void evolve(population &) const;
		std::string human_readable_extra() const;
	private:
		struct nlopt_wrapper_data
		{
			problem::base const		*prob;
			decision_vector			*x;
			fitness_vector			*f;
			constraint_vector		*c;
			problem::base::c_size_type	c_comp;
		};
		static double objfun_wrapper(int, const double *, double *, void *);
		static double constraints_wrapper(int, const double *, double *, void *);
	private:
		const nlopt_algorithm	m_algo;
		const bool		m_constrained;
		const std::size_t	m_max_iter;
		const double		m_tol;
};

}}

#endif
