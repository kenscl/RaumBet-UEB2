#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include "helper/mathhelper.h"
#include "sgp4/coordinates.h"
#include "sgp4/timeDate.h"
#include "tle/tle.h"
#include "tle/tlereader.h"
#include "sgp4/timeDate.h"
#include "sgp4/coordinates.h"
#include "sgp4/SGP4Propagator.h"
#include <iomanip>

#ifndef __MAIN 
#define __MAIN
int main(int argc, char *argv[])
{
	// File with TLEs
	string fileName = "../tle.txt";

	// Get Map from txt csv_file
	map<int, TLE> satMap = readTlesFromFile(fileName.c_str());

	TLE test_case_tle = satMap[88888];	
	TLE sonate = satMap[59112];
    ECICoordinate POS, VEL;
    SGP4Propagator propagator;

    // Generate csv-file from Sonate2 orbit for 100 min
    std::ofstream csv_file("sonate2_100min_orbit.csv");
    if (!csv_file.is_open()) {
        printf("Datei konnte nicht geöffnet werden!\n");
        return 1;
    }

    ECICoordinate pos, vel;
    GeocentricCoordinate eciToGeocentric;
    GeodeticCoordinate eciToGeodetic;
    propagator.setTle(sonate);
    int32_t secsAfterEpoch = 0;
    double jd, gmst = 0;
    int standardPrecision = 2, highPrecision = 8;       // amount of decimal places

    // loops in 5 min steps from 0 min to 100 min
    for (int i = 0; i < 100; i += 5) {
        secsAfterEpoch = 60 * i;
        propagator.calculatePositionAndVelocity(secsAfterEpoch, pos, vel);
        jd = computeJD(sonate.getYear(), sonate.getDayFraction() + secsAfterEpoch / 86400.0);
        gmst = computeGMST(jd);
        eciToGeocentric = convertECItoGeocentric(pos, jd);
        eciToGeodetic = convertECItoGeodetic(pos, jd);

        // time after epoch [min]
        csv_file << 5 * i << ",";
        // ECI coordinates of satellite [km]
        csv_file << std::fixed << std::setprecision(standardPrecision) << pos.x << "," << pos.y << "," << pos.z << ",";
        // current Julian date [days]
        csv_file << std::setprecision(highPrecision) << jd << std::setprecision(standardPrecision) << ",";
        // Greenwich Mean Sidereal Time [deg]
        csv_file << rad2deg(gmst) << ",";
        // Geocentric coordinates of satellite
        csv_file << rad2deg(eciToGeocentric.latitude) << "," << rad2deg(eciToGeocentric.longitude) << "," << eciToGeocentric.hight << ",";
        // Geodetic coordinates of satellite
        csv_file << rad2deg(eciToGeodetic.latitude) << "," << rad2deg(eciToGeodetic.longitude) << "," << eciToGeodetic.hight << "\n";
    }

    csv_file.close();

    /*printf("Test Case: \n");
    propagator.setTle(test_case_tle);
    for (int i = 0; i <= 1440 * 60; i += 360 * 60) {
        propagator.calculatePositionAndVelocity(i, POS, VEL);
        printf("Time after Epoch %d \n", i / 60);
        printf("Position: %f %f %f \n", POS.x, POS.y, POS.z);
        printf("Velocity: %f %f %f \n \n", VEL.x, VEL.y, VEL.z);

    }
    printf("--------------------------------------------------------------------\n");
    printf("Sonate: \n");
    propagator.setTle(sonate);
    for (int i = 0; i <= 1440 * 60; i += 360 * 60) {
        propagator.calculatePositionAndVelocity(i, POS, VEL);
        printf("Time after Epoch %d \n", i / 60);
        printf("Position: %f %f %f \n", POS.x, POS.y, POS.z);
        printf("Velocity: %f %f %f \n \n", VEL.x, VEL.y, VEL.z);
    }
    printf("--------------------------------------------------------------------\n");
	*/
    /*ECICoordinate eci;
    eci.x = 1000000;
    eci.y = 1000000;
    eci.z = 1000000;
    printf ("std: %f \n", computeJD(24,133.5));
    GeocentricCoordinate gd = convertECItoGeocentric(eci, computeJD(24, 133.5));
    Vector_3D ecef;
    ecef = geozentricToECEF(gd);
    convertECItoECEF(eci, computeJD(24,1.5)).print();
    ecef.print();
    printf ("gd %f %f %f \n", gd.latitude, gd.longitude, gd.hight);*/

    //propagator.setTle(sonate); for (int i = 0; i < 100 * 60; i++){ propagator.calculatePositionAndVelocity(i, POS, VEL);
    //    double epoch = computeJD(24, sonate.getDayFraction() + (double) i / (24 * 60 * 60));
    //    GeocentricCoordinate gc = convertECItoGeocentric(POS, epoch);
    //    GeodeticCoordinate gd = convertECItoGeodetic(POS, epoch);
    //    printf("%d, %f, %f, %f \n",i, rad2deg(gc.latitude - gd.latitude) , rad2deg(gc.longitude - gd.longitude), gc.hight - gd.hight);
    //}
    //

    propagator.setTle(sonate);
    printf("mm %f \n", sonate.getMeanMotion());
    for (int i = 0; i < 3 * 60 * (TWO_PI / sonate.getMeanMotion()); i++) {
        propagator.calculatePositionAndVelocity(i, POS, VEL);
        double epoch = computeJD(24, sonate.getDayFraction() + (double) i / (24 * 60 * 60));
        GeodeticCoordinate gd = convertECItoGeodetic(POS, epoch);
        printf("%f %f \n", rad2deg(gd.longitude), rad2deg(gd.latitude));
    }

}
#endif
