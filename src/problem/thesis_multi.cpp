#include <string>
#include <vector>
#include <numeric>
#include <cmath>

#include "thesis_multi.h"
#include "/home/enne/sml/include/sml.hpp"
#include "/home/enne/pykep/src/lambert_problem.h"


namespace pagmo { namespace problem {


/// Constructor
/**
 * Constructs the problem
 *
 * @param[in]
 * @param[in] 
 */
thesis_multi::thesis_multi( const unsigned int aDim, 
                            const unsigned int aNumberOfLegs,
                            const double aStayTime,
                            const double aDepartureEpochUpperBound,
                            const double aTimeOfFlightUpperBound,
                            const DateTime anInitialEpoch,
                            std::vector< Tle > aTleList)
                    :   base( aDim ), 
                        numberOfLegs( aNumberOfLegs ),
                        stayTime( aStayTime ),
                        departureEpochUpperBound( aDepartureEpochUpperBound ),
                        timeOfFlightUpperBound( aTimeOfFlightUpperBound ),
                        initialEpoch( anInitialEpoch ),
                        TleList( aTleList )
{
    // Construct here the problem (bounds etc.)
    std::vector<double>::size_type dim(3*numberOfLegs+1);
    std::vector<double> lb(dim);
    std::vector<double> ub(dim);

    for (unsigned int i = 0; i < numberOfLegs+1; ++i)
    {
        lb[i] = 0;
        ub[i] = static_cast< int >(TleList.size());
    }

    for (unsigned int i = 0; i < numberOfLegs; ++i)
    {
        lb[numberOfLegs+1+i*2] = 0.0;
        lb[numberOfLegs+1+1+i*2] = 1.0;
        ub[numberOfLegs+1+i*2] = departureEpochUpperBound;
        ub[numberOfLegs+1+1+i*2] = timeOfFlightUpperBound;
    }
    set_bounds(lb,ub);
}


/// Clone method.
base_ptr thesis_multi::clone() const
{
    return base_ptr(new thesis_multi(*this));
}

/// Implementation of the objective function.
void thesis_multi::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
    std::vector<int> integerSequence;
    std::vector< Tle > sequence;
    for (unsigned int k = 0; k < numberOfLegs+1; ++k)
    {
        int currentInteger = static_cast< int >(x[k]);
        integerSequence.push_back(currentInteger);
        sequence.push_back(TleList[currentInteger]);
    }

    std::sort( integerSequence.begin( ), integerSequence.end( ) );           
    if (adjacent_find(integerSequence.begin( ),integerSequence.end( ) ) != integerSequence.end() )
    {
        f[0]=55555555;
    }
    else
    {
        // Initialize vector to store dV of each leg
        std::vector<double> deltaVs(numberOfLegs,1000);
        // Set previous
        DateTime previousArrival = initialEpoch;

        for (unsigned int leg = 0; leg < numberOfLegs; ++leg)
        {
            double currentDepartureEpoch = x[numberOfLegs+1+2*leg];
            double currentTimeOfFlight = x[numberOfLegs+2+2*leg];
            Tle departureObject    = sequence[0+leg];
            Tle arrivalObject      = sequence[1+leg];
    
            DateTime departureEpoch = initialEpoch.AddSeconds(currentDepartureEpoch);
                
            if (previousArrival > departureEpoch)
            {
                deltaVs[leg] = 222222;
            }
            else
            {

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
                DateTime arrivalEpoch = departureEpoch.AddSeconds( currentTimeOfFlight );
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
                                                       currentTimeOfFlight, // timeOfFlight
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
                
                // Store lowest DV for outside for loop
                deltaVs[leg] = *minimumDeltaVIterator;             
                
                // Save arrival epoch plus stay time for check in next iteration
                previousArrival = arrivalEpoch.AddSeconds(stayTime);
            }
        }
        f[0] = std::accumulate(deltaVs.begin(), deltaVs.end(), 0.0);
    }
}


thesis_multi::~thesis_multi() {} 

/// Gets the data member
/**
 * @return[out] data member
 */
// const double& thesis_multi::get_member() const {
//     return m_member;
// }


/// Gets the departure object
/**
 * @return[out] departure object
 */
// const Tle& thesis_multi::get_departure_object() const {
//     return departureObject;
// }


/// Gets the arrival object
/**
 * @return[out] arrival object
 */
// const Tle& thesis_multi::get_arrival_object() const {
//     return arrivalObject;
// }

/// Sets the data member
/**
 *
 * @param[in] value new value
//  */
// void thesis_multi::set_member(double value) {
//     m_member = value;
// }

// std::string thesis_multi::a_method() const {
//     return "I do null";
// }


std::string thesis_multi::get_name() const {
    return "Thesis_multi Enne Joepi de poepie!!";
}


/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the values vector, the weights vectors and the max weight. It is concatenated
 * with the base::problem human_readable
//  */
// std::string thesis_multi::human_readable_extra() const
// {
//     std::ostringstream oss;
//     oss << "\n\tData member value: " << m_member;
//     return oss.str();
// }

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::thesis_multi)
