#ifndef PAGMO_PROBLEM_THESIS_H
#define PAGMO_PROBLEM_THESIS_H

#include <string>
#include "../config.h"
#include "../exceptions.h"
#include "../serialization.h"
#include "../types.h"
#include "../problem/base.h"
#include "../../../sgp4/libsgp4/Tle.h"
#include "../../../sgp4/libsgp4/SGP4.h"

// #include <libsgp4/Eci.h>
// #include <libsgp4/Globals.h>
// #include <libsgp4/Tle.h>
// #include <libsgp4/SGP4.h>


// #include <pagmo/src/config.h>
// #include <pagmo/src/serialization.h>
// #include <pagmo/src/types.h>
// #include <pagmo/src/problem/base.h>

namespace pagmo{ namespace problem {

/// A random problem
/**
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class __PAGMO_VISIBLE thesis: public base
{
    public:
        thesis( unsigned int aDim       = 2, 
                double departureEpochUpperBound = 86400,
                double timeOfFlightUpperBound = 86400,
                const Tle departureObject    = Tle(
                    "1 22830U 93061H   16010.66911318 -.00000005  00000-0  15512-4 0  9991",
                    "2 22830  98.8796   7.9979 0011279  21.4730 107.6291 14.31392246163351"),
                const Tle arrivalObject      = Tle(
                    "1 20443U 90005H   16011.52874183  .00000186  00000-0  77220-4 0  9991", 
                    "2 20443  98.7809  34.7119 0011411 108.4586  53.3704 14.37782899359716"),
                const DateTime initialEpoch  =  DateTime( 2016,1,12,12,0,0 )
                );

        base_ptr clone() const;
        std::string get_name() const;

        std::string a_method() const;
        void set_member(const double);

    protected:
        void objfun_impl(fitness_vector &, const decision_vector &) const;
        // std::string human_readable_extra() const;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base>(*this);
        }
        double departureEpochUpperBound;
        double timeOfFlightUpperBound;
        const Tle departureObject;
        const Tle arrivalObject;
        DateTime initialEpoch;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::thesis)
#endif // PAGMO_PROBLEM_THESIS_H
