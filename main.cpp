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

	//TLE tle = satMap[12345];	// Beispiel


	// Aufgabe: Ergänzen Sie hier ihren Code


}
#endif
