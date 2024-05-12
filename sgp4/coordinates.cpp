#include "coordinates.h"

Vector_3D convertECItoECEF(const ECICoordinate &eciCoord, double jd){
    Matrix_3D rot;
    Vector_3D eci, ecef;

    eci.x = eciCoord.x;
    eci.y = eciCoord.y;
    eci.z = eciCoord.z;
    
    rot.r[0][0] = cos (jd);
    rot.r[0][1] = sin (jd);
    rot.r[0][2] = 0;

    rot.r[1][0] = -sin (jd);
    rot.r[1][1] = cos (jd);
    rot.r[1][2] = 0;

    rot.r[2][0] = 0;
    rot.r[2][1] = 0;
    rot.r[2][2] = 1;

    ecef = rot * eci;
    return ecef;

}

GeocentricCoordinate convertECItoGeocentric(const ECICoordinate &eciCoord, double jd){
    Vector_3D ecef;
    ecef = convertECItoECEF(eciCoord, jd);

    GeocentricCoordinate result;
    result.hight = sqrt (ecef.x * ecef.x + ecef.y * ecef.y + ecef.z * ecef.z);
    result.longitude= atan2 (ecef.z, sqrt (ecef.x * ecef.x + ecef.y * ecef.y));
    result.latitude = atan2 (ecef.y, ecef.x);
    return result;
}

GeodeticCoordinate convertECItoGeodetic(const ECICoordinate &eciCoord, double jd){
    GeocentricCoordinate geocentric;
    GeodeticCoordinate geodetic;
    geodetic.latitude = geocentric.latitude;
    geodetic.hight = geocentric.hight;
    geocentric = convertECItoGeocentric(eciCoord, jd);
    geodetic.longitude= 1/298.257223563 * tan (geocentric.longitude);

    return geodetic;
}
