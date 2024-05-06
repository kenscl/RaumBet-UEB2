#ifndef COORDINATES_H_
#define COORDINATES_H_


struct GeocentricCoordinate {
    double x,y,z;
};

struct GeodeticCoordinate {
    double x,y,z;
};

/**
 * @brief Simple 3D Cartesian coordinate in ECI
 */
struct ECICoordinate {
    double x,y,z;
};


/**
 * @brief Converts the ECI coordinate to the geodetic system
 */
GeodeticCoordinate convertECItoGeodetic(const ECICoordinate &eciCoord, double jd);

/**
 * @brief Converts the ECI coordinate to the geocentric system
 */
GeocentricCoordinate convertECItoGeocentric(const ECICoordinate &eciCoord, double jd);

#endif /* COORDINATES_H_ */
