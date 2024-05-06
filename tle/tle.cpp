#define _USE_MATH_DEFINES
#include "tle.h"
#include "../helper/mathhelper.h"
#include <iostream>
#include <iomanip>


/* -------------------------- Constructor --------------------- */

TLE::TLE() : valid(false)
{

}

/**
 * @brief Creates a new two line elements object from the raw data contained in the two lines
 *
 * First checks if the checksum is valid and line1 and line2 have the correct length.
 * Then parses the two lines into the class variables.
 *
 * @param line0	the first line of the TLE (contains the name of the object)
 * @param line1 the second line of the TLE (contains the deviation of the orbital elements + epoch)
 * @param line2 the third line of the TLE (contains the orbital elements)
 */
TLE::TLE(const char *line0, const char *line1, const char *line2) : TLE()
{
    valid = populate(line0, line1, line2);
}

/* --------------------------- public functions ------------------------ */

/**
 * @brief Prints the TLE variables to console
 * 
 * Prints variables if TLE is valid, an 'invalid' notifier otherwise.
 */
void TLE::print()
{
    if (valid) {
		printf("Name:             %s\n", satelliteName);
        printf("SatNumber:        %d\n", satelliteNr);
        printf("Designator:       %s\n", intDesignator);
        printf("Year:             %d\n", year);
        printf("DayFrac:          %10.6f\n", dayFraction);
        printf("bStar:            %10.6f\n", bStar);
        printf("Inclination:      %10.6f [rad]\t%6.4f [deg]\n", inclination, rad2deg(inclination));
        printf("RAAN:             %10.6f [rad]\t%6.4f [deg]\n", raan, rad2deg(raan));
        printf("Eccentricity:     %10.6f\n", eccentricity);
        printf("ArgOfPerigee:     %10.6f [rad]\t%6.4f [deg]\n", argumentOfPerigee, rad2deg(argumentOfPerigee));
        printf("MeanAnomaly:      %10.6f [rad]\t%6.4f [deg]\n", meanAnomaly, rad2deg(meanAnomaly));
        printf("MeanMotion:       %10.6f [rad/min]\t%6.4f [deg/min]\n", meanMotion, rad2deg(meanMotion));
		printf("------------------------------------------------------\n");
		printf("True Anomaly:     %10.6f [rad]\t%6.4f [deg]\n", calcTrueAnomaly(), rad2deg(calcTrueAnomaly()));
		printf("Semi-major axis: %10.6f [km]\n", calcSemiMajorAxis());
		printf("------------------------------------------------------\n");

	} else {
        printf(" TLE INVALID!\n");
    }
}

/**
* @brief Calculates and returns the orbit period
* @return orbit period [s]
*/
double TLE::calcOrbitPeriod()
{
	double period = 2 * M_PI / meanMotion;                  // [min]
	period *= 60;                                           // [s]
	return period;
}

/**
 * @brief Calculates and returns the semi-major axis
 * 
 * First, calculates the period from the mean motion. As the mean motion is in [rad/min],
 * the period has to be converted to [s].
 * Then, the semi-major axis is calculated with the period and mu_Earth.
 *
 * @return semi-major axis
 */
double TLE::calcSemiMajorAxis()
{
    // Zur Berechnung siehe auch RFBA_2.pdf, Seite 15
    double mu = 398600;                                     // [km³/s²]
    double period = 2 * M_PI / meanMotion;                  // [min]
    period *= 60;                                           // [s]
    double sma = cbrt(mu *period*period / (4*M_PI*M_PI) );  // [km]

    return sma;
}

/**
 * @brief Calculates and returns the true anomaly
 *
 * Uses Newton's method to determine the eccentric anomaly. Starting value is the mean anomaly.
 * Four iterations are usually sufficient.
 * The true anomaly is calculated from Kepler's equation. The result is normalized to [0,2*pi]
 *
 * @return true anomaly
 */
double TLE::calcTrueAnomaly()
{
    // Zur Berechnung siehe auch RFBA_2.pdf, Seite 15
    // Eccentric Anomaly

	std::cout << std::setprecision(18);

    double E = meanAnomaly;

    double g, gdot;
    for (int i = 0; i< 4; i++) {
        g = E - eccentricity * sin(E) - meanAnomaly;
        gdot = 1 - eccentricity * cos(E);
        E = E - (g/gdot);

		/*std::cout << "g = \t" << g << std::endl;
		std::cout << "gdot = \t" << gdot << std::endl;
		std::cout << "E = \t" << E << std::endl;*/

    }

    // True Anomaly
    double trueAnomaly = 2 * atan2( tan(E/2) * sqrt( (1+eccentricity) / (1-eccentricity) ) , 1);
	return (trueAnomaly < 0 ? trueAnomaly + TWO_PI : trueAnomaly);
}

/* ----------------------- private functions ---------------------------- */

/**
 * @brief Checks if a TLE line is valid.
 *
 * Therefore checks if the TLE has exactly 69 characters and the
 * checksum of the TLE is correct.
 *
 * @param line	the TLE line to check
 *
 * @return true if the TLE is valid, otherwise false
 */
bool TLE::isLineValid(const char *line) const
{
	//check if the TLE has exactly 69 characters
	int32_t len = strlen(line);
	if (len != 69) {
		return false;
	}

	//check if the checksum is valid
	return isChecksumValid(line);
}

/**
 * @brief Checks if the checksum of the TLE line is correct
 *
 * @param line	the TLE line to check
 *
 * @return true if the checksum is correct, otherwise false
 */
bool TLE::isChecksumValid(const char *line) const
{
	//compute the checksum
	int32_t cs = 0;
	for (int i = 0; i < 68; ++i) {
		if (isdigit(line[i])) {
			cs += (line[i] - '0');
		} else if (line[i] == '-') {
			cs += 1;
		}
	}

    //modulo 10 of checksum
	cs %= 10;

	return (line[68] - '0') == cs;
}

/**
 * @brief Populates the class variables from a TLE.
 *
 * First checks if the TLE is valid. Then parses the TLE into the
 * class variables.
 *
 * @param line0	the first line of the TLE (contains the name of the object)
 * @param line1 the second line of the TLE (contains the deviation of the orbital elements + epoch)
 * @param line2 the third line of the TLE (contains the orbital elements)
 *
 * @return true if the TLE was parsed, false otherwise
 */
bool TLE::populate(const char *line0, const char *line1, const char *line2)
{
	//check if the TLE is valid before parsing it
	bool tleValid = isLineValid(line1);
	tleValid &= isLineValid(line2);

	//return if the TLE is invalid
	if (!tleValid) {
		return false;
	}

    //TLE is valid beyond this point
    strncpy(satelliteName, line0, 24);
    for (int i = 23; i >= 0; --i) {
        if (isblank(satelliteName[i])){
            satelliteName[i] = '\0';
        } else {
            break;
        }
    }

	//copy the satellite number
    satelliteNr = strtol(&line1[2], nullptr, 10);

	//copy the international designator
	strncpy(intDesignator, &line1[9], 8);
    for (int i = 7; i>0; --i) {
        if (!isblank(intDesignator[i])){
            break;
        } else {
           intDesignator[i] = '\0';
        }
    }

	//find year of tle
	char temp[12] = { '\0' };
	strncpy(temp, &line1[18], 2);
	year = strtol(temp, nullptr, 10);
	if (year <= 56) {
		year += 2000;
	} else {
		year += 1900;
	}

    char *pEnd;
	//parse fraction of the year
    dayFraction = strtod(&line1[20], &pEnd);

    if (*pEnd == '.') {
        char *pEndNew = nullptr;
        double afterComma = strtod(&pEnd[1], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 1);
        for (int i = 0; i < exponent; ++i) {
            afterComma *= 0.1;
        }
        dayFraction += afterComma;
    }

	//parse first derivative of mean motion
    nDot = strtod(&line1[33], &pEnd);
    if (*pEnd == '-') {
        char *pEndNew = nullptr;
        nDot = -1;
        nDot *= strtod(&line1[35], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 2);
        for (int i = 0; i < exponent; ++i) {
            nDot *= 0.1;
        }

    } else if (*pEnd == ' ' || *pEnd == '+') {
        char *pEndNew = nullptr;
        nDot = 1.0;
        nDot *= strtod(&line1[35], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 2);
        for (int i = 0; i < exponent; ++i) {
            nDot *= 0.1;
        }
    }

	//parse second derivative of mean motion
	strncpy(temp, &line1[44], 6);
	int32_t tempInt = strtol(temp, nullptr, 10);
	int32_t exponent = strtol(&line1[51], nullptr, 10);
	nDotDot = tempInt * 1.0e-5;
	for (int i = 0; i < exponent; ++i) {
		nDotDot *= 0.1;
	}


	//parse Bstar drag term
	strncpy(temp, &line1[53], 6);
	tempInt = strtol(temp, nullptr, 10);
	exponent = strtol(&line1[60], nullptr, 10);
	bStar = tempInt * 1.0e-5;
	for (int i = 0; i < exponent; ++i) {
		bStar *= 0.1;
    }

	//second line of the TLE
	//parse orbital elements
    inclination = strtod(&line2[8], &pEnd);
    if (*pEnd == '.') {
        char *pEndNew = nullptr;
        double afterComma = strtod(&pEnd[1], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 1);
        for (int i = 0; i < exponent; ++i) {
            afterComma *= 0.1;
        }
        inclination += afterComma;
    }
    inclination = deg2rad(inclination);

    raan = strtod(&line2[17], &pEnd);

    if (*pEnd == '.') {
        char *pEndNew = nullptr;
        double afterComma = strtod(&pEnd[1], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 1);
        for (int i = 0; i < exponent; ++i) {
            afterComma *= 0.1;
        }
        raan += afterComma;
    }
	raan = deg2rad(raan);

    eccentricity = strtol(&line2[26], nullptr, 10);
	eccentricity *= 0.0000001;

    argumentOfPerigee = strtod(&line2[34], &pEnd);
    if (*pEnd == '.') {
        char *pEndNew = nullptr;
        double afterComma = strtod(&pEnd[1], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 1);
        for (int i = 0; i < exponent; ++i) {
            afterComma *= 0.1;
        }
        argumentOfPerigee += afterComma;
    }
	argumentOfPerigee = deg2rad(argumentOfPerigee);

    meanAnomaly = strtod(&line2[43], &pEnd);
    if (*pEnd == '.') {
        char *pEndNew = nullptr;
        double afterComma = strtod(&pEnd[1], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 1);
        for (int i = 0; i < exponent; ++i) {
            afterComma *= 0.1;
        }
        meanAnomaly += afterComma;
    }
	meanAnomaly = deg2rad(meanAnomaly);

	//copy mean motion
	strncpy(temp, &line2[51], 11);
    meanMotion = strtod(temp, &pEnd);
    if (*pEnd == '.') {
        char *pEndNew = nullptr;
        double afterComma = strtod(&pEnd[1], &pEndNew);
        int32_t exponent = (pEndNew - pEnd - 1);
        for (int i = 0; i < exponent; ++i) {
            afterComma *= 0.1;
        }
        meanMotion += afterComma;
    }

    meanMotion *= (TWO_PI / 1440.0);

	//copy revolutions at epoch
	strncpy(temp, &line2[63], 5);
	temp[5] = '\0';
	revsAtEpoch = strtol(temp, nullptr, 10);

	return true;
}
