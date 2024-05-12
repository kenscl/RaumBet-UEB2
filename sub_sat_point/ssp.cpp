#include "ssp.h"
#include <cmath>

void sub_satellite_point_round_earth(ECICoordinate location, GeocentricCoordinate &on_earth){
    on_earth.breite = atan2(location.z, sqrt(location.x * location.x + location.y * location.y));
}
