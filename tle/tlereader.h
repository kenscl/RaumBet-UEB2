#ifndef TLEREADER_H_
#define TLEREADER_H_

#include <map>
#include "tle.h"

using namespace std;

/**
 * @brief Reads TLEs from a txt file and saves them in a map with the satellite number as key.
 * @param fileName - txt file with TLEs
 * @return Map of Satellite numbers to their TLEs
 */
map<int, TLE> readTlesFromFile(const char *fileName);

#endif /* TLEREADER_H_ */