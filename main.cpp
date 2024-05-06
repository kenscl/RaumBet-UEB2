#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include "tle/tle.h"
#include "tle/tlereader.h"
#include "sgp4/timeDate.h"
#include "sgp4/coordinates.h"
#include "sgp4/SGP4Propagator.h"

#ifndef __MAIN 
#define __MAIN
int main(int argc, char *argv[])
{
	// File with TLEs
	string fileName = "../tle.txt";

	// Get Map from txt file
	map<int, TLE> satMap = readTlesFromFile(fileName.c_str());

	TLE test_case_tle = satMap[88888];	
	TLE sonate = satMap[59112];	
    ECICoordinate POS;
    ECICoordinate VEL;
    SGP4Propagator propagator;
    
    printf("Test Case: \n");
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


}
#endif
