#ifndef _sgp4_sgp4Propagator_
#define _sgp4_sgp4Propagator_

#include "../tle/tle.h"
#include "coordinates.h"


/**
 * @brief This class implements the SGP4 orbit propagator
 *
 * @author Andreas Maurer, andreas.maurer@uni-wuerzburg.de
 */
class SGP4Propagator
{
private:
	TLE tle; //< contains the TLE used for the orbit propagation
	bool neConstsInitialized = false; ///< indicates if the near earth constants were initialized

public:
	SGP4Propagator();

	/**
	 * @brief Sets the TLE used for the orbit propagation with the SGP4 model
	 *
	 * This function is also used to update the TLE once a new one is available.
	 */
	inline void setTle(const TLE &_tle)
	{
		//set the TLE to the new value
		this->tle = _tle;

		//mark initialization of near earth constant as NOT done
		neConstsInitialized = false;
	}

	/**
	 * @return true if the NearEarthConstant are initialized
	 */
	inline bool areNeConstsInitialized() const
	{
		return neConstsInitialized;
	}

	/**
	 * @brief Calculates the position and velocity of the satellite in the ECI coordinate frame using the SGP4 model
	 *
	 * @param secsAfterEpoch	Seconds after the epoch of the TLE [s]
	 * @param satPos			Position of the satellite in ECI [km]
	 * @param satVel			Velocity of the satellite in ECI [km/s]
	 */
	void calculatePositionAndVelocity(int32_t secsAfterEpoch, ECICoordinate &satPos, ECICoordinate &satVel);
};

#endif /* end _sgp4_sgp4Propagator_ */

