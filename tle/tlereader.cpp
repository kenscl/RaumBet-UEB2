#include "tlereader.h"
#include <iostream>
#include <fstream>
#include <string>

map<int, TLE> readTlesFromFile(const char *fileName) {

	map<int, TLE> satMap;

	// Open file
	ifstream tleFile;
	tleFile.open(fileName);

	// Strings to hold TLE lines
	string line0;
	string line1;
	string line2;

	int cnt_valid = 0, cnt_invalid = 0;

	cout << "Reading TLEs from " << fileName << endl << endl;

	// Iterate over txt file until EOF
	while (tleFile.good()) {
		// Only if all three lines found
		if (getline(tleFile, line0)
			&& getline(tleFile, line1)
			&& getline(tleFile, line2)) {

			// Catches problem on linux systems: Tle txt files contain windows line breaks. '\r' is left over from "\r\n".
			if (line0.back() == '\r') line0.pop_back();
			if (line1.back() == '\r') line1.pop_back();
			if (line2.back() == '\r') line2.pop_back();

			// Construct TLE from lines
			TLE tle(line0.data(), line1.data(), line2.data());

			if (!tle.isValid()) {
				cout << "Skipped Invalid TLE!" << endl;
				cnt_invalid++;
			}
			else {
				// Add to map
				satMap[tle.getSatelliteNr()] = tle;
				//cout << tle.getSatelliteName() << " added to Map." << endl;
				cnt_valid++;
			}
		}
	}
	cout << endl << "EOF reached. Finished." << endl;
	cout << "Added " << cnt_valid << " valid, " << cnt_invalid << " invalid TLEs." << endl << endl;

	tleFile.close();
	return satMap;
}
