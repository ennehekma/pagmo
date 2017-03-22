#include <string>
#include <vector>
#include <numeric>
#include <cmath>

#include "thesis.h"
#include <SML/sml.hpp>
#include <keplerian_toolbox/lambert_problem.h>

namespace pagmo { namespace problem {


/// Constructor
/**
 * Constructs the problem
 *
 * @param[in] seq std::vector of kep_toolbox::planet_ptr containing the encounter sequence for the trajectoty (including the initial planet)
 * @param[in] t0_l kep_toolbox::epoch representing the lower bound for the launch epoch
 */
thesis::thesis( unsigned int aDim, 
                // double aData)
                // double aData,
                // double hoeDan,
                const Tle aDepartureObject,
                const Tle anArrivalObject,
                const DateTime anInitialEpoch,
                const double aDepartureEpochUpperBound,
                const double aTimeOfFlightUpperBound)
    :   base(aDim), 
        // m_member(aData),
        // m_boe(hoeDan),
        departureObject( aDepartureObject ),
        arrivalObject( anArrivalObject ),
        initialEpoch (anInitialEpoch),
        departureEpochUpperBound (aDepartureEpochUpperBound),
        timeOfFlightUpperBound (aTimeOfFlightUpperBound)

{
    // Construct here the problem (bounds etc.)
    std::vector<double>::size_type dim(2);
    std::vector<double> lb(dim);
    std::vector<double> ub(dim);

    lb[0] = 0;
    lb[1] = 1;
    // lb[2] = 0;
    // lb[3] = 0;
    // lb[4] = 0;
    
    ub[0] = departureEpochUpperBound;
    ub[1] = 2*86400;
    // ub[2] = 100;
    // ub[3] = 100;
    // ub[4] =   1;
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
    // Define departure position:
    const DateTime departureEpoch = initialEpoch.AddSeconds(x[0]);
    SGP4 sgp4Departure( departureObject );
    const Eci tleDepartureState = sgp4Departure.FindPosition( departureEpoch );

    boost::array< double, 3 > departurePosition;
    departurePosition[ 0 ] = tleDepartureState.Position( ).x;
    departurePosition[ 1 ] = tleDepartureState.Position( ).y;
    departurePosition[ 2 ] = tleDepartureState.Position( ).z;
    
    boost::array< double, 3 > departureVelocity;
    departureVelocity[ 0 ] = tleDepartureState.Velocity( ).x;
    departureVelocity[ 1 ] = tleDepartureState.Velocity( ).y;
    departureVelocity[ 2 ] = tleDepartureState.Velocity( ).z;

    // Define arrival position:
    const DateTime arrivalEpoch = departureEpoch.AddSeconds( x[1] );
    SGP4 sgp4Arrival( arrivalObject );
    const Eci tleArrivalState   = sgp4Arrival.FindPosition( arrivalEpoch );

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

    const int numberOfSolutions = targeter.get_v1( ).size( );

    // Compute Delta-Vs for transfer and determine index of lowest.
    typedef std::vector< boost::array< double, 3 > > VelocityList;
    VelocityList departureDeltaVs( numberOfSolutions );
    VelocityList arrivalDeltaVs( numberOfSolutions );

    typedef std::vector< double > TransferDeltaVList;
    TransferDeltaVList transferDeltaVs( numberOfSolutions );

    for ( int i = 0; i < numberOfSolutions; i++ )
    {
        // Compute Delta-V for transfer.
        const boost::array< double, 3 > transferDepartureVelocity = targeter.get_v1( )[ i ];
        const boost::array< double, 3 > transferArrivalVelocity = targeter.get_v2( )[ i ];

        departureDeltaVs[ i ] = sml::add( transferDepartureVelocity,
                                          sml::multiply( departureVelocity, -1.0 ) );
        arrivalDeltaVs[ i ]   = sml::add( arrivalVelocity,
                                          sml::multiply( transferArrivalVelocity, -1.0 ) );

        transferDeltaVs[ i ]
            = sml::norm< double >( departureDeltaVs[ i ] )
                + sml::norm< double >( arrivalDeltaVs[ i ] );
    }
    const TransferDeltaVList::iterator minimumDeltaVIterator
        = std::min_element( transferDeltaVs.begin( ), transferDeltaVs.end( ) );
    // const int minimumDeltaVIndex
    //     = std::distance( transferDeltaVs.begin( ), minimumDeltaVIterator );

    f[0] = *minimumDeltaVIterator;
    // std::cout << f[0] << std::endl;
    
    // f[0] = 0.0;
    // for (decision_vector::size_type i = 0; i < x.size(); ++i) {
    //     f[0] += x[i];
    // }
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
