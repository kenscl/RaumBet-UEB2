/*
 * tle.h
 *
 *  Created on: 17.07.2019
 *      Author: MaurerAndreas
 */

#ifndef TLE_H_
#define TLE_H_

#include <stdint.h>
#include <string>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <ctype.h>
#include <cmath>
#include "../helper/mathhelper.h"

class TLE
{
private:
    char satelliteName[25] = {'\0'};
    int32_t satelliteNr;
    char intDesignator[9] = {'\0'};
	int32_t year;
	double dayFraction;
	double nDot;
	double nDotDot;
	double bStar;

	double inclination; ///< inclination [rad]
	double raan; ///< right ascension of ascending node [rad]
	double eccentricity; ///< eccentricity of the orbit
	double argumentOfPerigee; ///< argument of perigee [rad]
	double meanAnomaly; ///< mean anomaly [rad]
	double meanMotion; ///< mean motion [rad/min]
	int32_t revsAtEpoch;

	bool valid = false;

	bool isLineValid(const char *line) const;
	bool isChecksumValid(const char *line) const;
	bool populate(const char *line0, const char *line1, const char *line2);

public:
	TLE();
	TLE(const char *line0, const char *line1, const char *line2);

    /**
     * @brief Prints TLE variables to console if TLE is valid. Otherwise prints "INVALID TLE".
     */
    void print();

    /**
     * @brief Calculates and returns the semi-major axis
     * @return semi-major axis [km]
     */
    double calcSemiMajorAxis();

    /**
     * @brief Calculates and returns the true anomaly
     * @return true anomaly [rad]
     */
    double calcTrueAnomaly();

	/**
	* @brief Calculates and returns the orbit period
	* @return orbit period [s]
	*/
	double calcOrbitPeriod();

	/**
	 * @return true if this TLE is valid, otherwise false
	 */
	inline bool isValid() const
	{
		return valid;
	}

    inline const char*  getSatelliteName() const
	{
		return satelliteName;
	}

    inline int32_t getSatelliteNr() const
	{
		return satelliteNr;
	}

    inline const char*  getIntDesignator() const {
		return intDesignator;
	}

	inline int32_t getYear() const
	{
		return year;
	}

	inline double getDayFraction() const
	{
		return dayFraction;
	}

	inline double getnDot() const
	{
		return nDot;
	}

	inline double getnDotDot() const
	{
		return nDotDot;
	}

	inline double getBstar() const
	{
		return bStar;
	}

	inline double getInclination() const
	{
		return inclination;
	}

	inline double getRaan() const
	{
		return raan;
	}

	inline double getEccentricity() const
	{
		return eccentricity;
	}

	inline double getArgumentOfPerigee() const
	{
		return argumentOfPerigee;
	}

	inline double getMeanAnomaly() const
	{
		return meanAnomaly;
	}

	inline double getMeanMotion() const
	{
		return meanMotion;
	}

	inline int32_t getRevolutionAtEpoch() const
	{
		return revsAtEpoch;
	}
};

#endif /* TLE_H_ */
