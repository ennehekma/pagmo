#include <string>
#include <vector>
#include <numeric>
#include <cmath>

#include "thesis.h"
#include <SML/sml.hpp>
#include <keplerian_toolbox/lambert_problem.h>
// #include "/home/enne/sml/include/sml.hpp"
// #include "/home/enne/pykep/src/lambert_problem.h"


namespace pagmo { namespace problem {


/// Constructor
/**
 * Constructs the problem
 *
 * @param[in] seq std::vector of kep_toolbox::planet_ptr containing the encounter sequence for the trajectoty (including the initial planet)
 * @param[in] t0_l kep_toolbox::epoch representing the lower bound for the launch epoch
 */
thesis::thesis( unsigned int aDim, 
                double aDepartureEpochUpperBound,
                double aTimeOfFlightUpperBound)
                const Tle aDepartureObject,
                const Tle anArrivalObject,
                const DateTime anInitialEpoch,
    :   base(aDim), 
        departureEpochUpperBound (aDepartureEpochUpperBound),
        timeOfFlightUpperBound (aTimeOfFlightUpperBound),
        departureObject( aDepartureObject ),
        arrivalObject( anArrivalObject ),
        initialEpoch (anInitialEpoch)

{
    // Construct here the problem (bounds etc.)
    std::vector<double>::size_type dim(2);
    std::vector<double> lb(dim);
    std::vector<double> ub(dim);

    lb[0] = 0;
    lb[1] = 1;
    
    ub[0] = departureEpochUpperBound;
    ub[1] = timeOfFlightUpperBound;
    set_bounds(lb,ub);
}


/// Clone method.
base_ptr thesis::clone() const
{
    return base_ptr(new thesis(*this));
}

/// Implementation of the objective function.
void thesis::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    
    DateTime departureEpoch = initialEpoch.AddSeconds(x[0]);
    SGP4 sgp4Departure( departureObject );
    Eci tleDepartureState = sgp4Departure.FindPosition( departureEpoch );

    boost::array< double, 3 > departurePosition;
    departurePosition[ 0 ] = tleDepartureState.Position( ).x;
    departurePosition[ 1 ] = tleDepartureState.Position( ).y;
    departurePosition[ 2 ] = tleDepartureState.Position( ).z;
    
    boost::array< double, 3 > departureVelocity;
    departureVelocity[ 0 ] = tleDepartureState.Velocity( ).x;
    departureVelocity[ 1 ] = tleDepartureState.Velocity( ).y;
    departureVelocity[ 2 ] = tleDepartureState.Velocity( ).z;

    // Define arrival position:
    DateTime arrivalEpoch = departureEpoch.AddSeconds( x[1] );
    SGP4 sgp4Arrival( arrivalObject );
    Eci tleArrivalState   = sgp4Arrival.FindPosition( arrivalEpoch );

    boost::array< double, 3 > arrivalPosition;
    arrivalPosition[ 0 ] = tleArrivalState.Position( ).x;
    arrivalPosition[ 1 ] = tleArrivalState.Position( ).y;
    arrivalPosition[ 2 ] = tleArrivalState.Position( ).z;

    boost::array< double, 3 > arrivalVelocity;
    arrivalVelocity[ 0 ] = tleArrivalState.Velocity( ).x;
    arrivalVelocity[ 1 ] = tleArrivalState.Velocity( ).y;
    arrivalVelocity[ 2 ] = tleArrivalState.Velocity( ).z;

    // Calculate minimum DV
    kep_toolbox::lambert_problem targeter( departurePosition,
                                           arrivalPosition,
                                           x[1], // timeOfFlight
                                           kMU, // earthGravitationalParameter,
                                           1,   // !input.isPrograde,
                                           50);// input.revolutionsMaximum 

    int numberOfSolutions = targeter.get_v1( ).size( );

    // Compute Delta-Vs for transfer and determine index of lowest.
    typedef std::vector< boost::array< double, 3 > > VelocityList;
    VelocityList departureDeltaVs( numberOfSolutions );
    VelocityList arrivalDeltaVs( numberOfSolutions );

    typedef std::vector< double > TransferDeltaVList;
    TransferDeltaVList transferDeltaVs( numberOfSolutions );

    for ( int i = 0; i < numberOfSolutions; i++ )
    {
        // Compute Delta-V for transfer.
        boost::array< double, 3 > transferDepartureVelocity = targeter.get_v1( )[ i ];
        boost::array< double, 3 > transferArrivalVelocity = targeter.get_v2( )[ i ];

        departureDeltaVs[ i ] = sml::add( transferDepartureVelocity,
                                          sml::multiply( departureVelocity, -1.0 ) );
        arrivalDeltaVs[ i ]   = sml::add( arrivalVelocity,
                                          sml::multiply( transferArrivalVelocity, -1.0 ) );

        transferDeltaVs[ i ]
            = sml::norm< double >( departureDeltaVs[ i ] )
                + sml::norm< double >( arrivalDeltaVs[ i ] );
    }
    TransferDeltaVList::iterator minimumDeltaVIterator
        = std::min_element( transferDeltaVs.begin( ), transferDeltaVs.end( ) );
    f[0] = *minimumDeltaVIterator;
}

/// Gets the data member
/**
 * @return[out] data member
 */
// const double& thesis::get_member() const {
//     return m_member;
// }


/// Gets the departure object
/**
 * @return[out] departure object
 */
// const Tle& thesis::get_departure_object() const {
//     return departureObject;
// }


/// Gets the arrival object
/**
 * @return[out] arrival object
 */
// const Tle& thesis::get_arrival_object() const {
//     return arrivalObject;
// }

/// Sets the data member
/**
 *
 * @param[in] value new value
//  */
// void thesis::set_member(double value) {
//     m_member = value;
// }

// std::string thesis::a_method() const {
//     return "I do null";
// }


std::string thesis::get_name() const {
    return "Thesis Enne Joepi de poepie!!";
}


/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight. It is concatenated
 * with the base::problem human_readable
//  */
// std::string thesis::human_readable_extra() const
// {
//     std::ostringstream oss;
//     oss << "\n\tData member value: " << m_member;
//     return oss.str();
// }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::thesis)
