#ifndef COORDINATES_H_
#define COORDINATES_H_
#include "../math/matlib.h"


struct GeocentricCoordinate {
    double hight, latitude, longitude;
};

struct GeodeticCoordinate {
    double hight, latitude, longitude;
};

/**
 * @brief Simple 3D Cartesian coordinate in ECI
 */
struct ECICoordinate {
    double x,y,z;
};

struct ECEFCoordinate {
    double x,y,z;
};


Vector_3D convertECItoECEF(const ECICoordinate &eciCoord, double jd);
/**
 * @brief Converts the ECI coordinate to the geodetic system
 */
GeodeticCoordinate convertECItoGeodetic(const ECICoordinate &eciCoord, double jd);

/**
 * @brief Converts the ECI coordinate to the geocentric system
 */
GeocentricCoordinate convertECItoGeocentric(const ECICoordinate &eciCoord, double jd);

#endif /* COORDINATES_H_ */
