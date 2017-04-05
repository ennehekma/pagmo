#ifndef PAGMO_PROBLEM_THESIS_MULTI_H
#define PAGMO_PROBLEM_THESIS_MULTI_H

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
class __PAGMO_VISIBLE thesis_multi: public base
{
    public:
        thesis_multi( 
            const unsigned int  aDim                        = 2, 
            const unsigned int  aNumberOfLegs               = 2,
            const double        stayTime                    = 0.0,
            const double        departureEpochUpperBound    = 86400.0,
            const double        timeOfFlightUpperBound      = 86400.0,
            const DateTime      initialEpoch                = DateTime( 2016,1,12,12,0,0 ),
            std::vector< Tle >  TleList                     = { Tle(
                "1 22830U 93061H   16010.66911318 -.00000005  00000-0  15512-4 0  9991",
                "2 22830  98.8796   7.9979 0011279  21.4730 107.6291 14.31392246163351"),
                Tle(
                "1 20443U 90005H   16011.52874183  .00000186  00000-0  77220-4 0  9991", 
                "2 20443  98.7809  34.7119 0011411 108.4586  53.3704 14.37782899359716")}
                );

        base_ptr clone() const;
        std::string get_name() const;

        std::string a_method() const;
        void set_member(const double);
	virtual ~thesis_multi();
        // const double& get_member() const;

    protected:
        void objfun_impl(fitness_vector &, const decision_vector &) const;
        // std::string human_readable_extra() const;

    private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int)
        {
            ar & boost::serialization::base_object<base>(*this);
            // ar & departureEpochUpperBound;
            // ar & timeOfFlightUpperBound;
            // ar & initialEpoch;
        }
        
        const unsigned int numberOfLegs;
        const double stayTime;
        const double departureEpochUpperBound;
        const double timeOfFlightUpperBound;
        const DateTime initialEpoch;
        std::vector< Tle > TleList;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::thesis_multi)
#endif // PAGMO_PROBLEM_THESIS_MULTI_H
