#include "timeDate.h"
#include <cmath>
#include <cstdio>
#include "../helper/mathhelper.h"

/*
 * Implementation of the algorithm presented in https://de.wikipedia.org/wiki/Julianisches_Datum.
 */
double computeJD(int year, double dayFraction) {
    double B, JD = 0;
    int month = 1;

    // Array to store days in each month, considering leap years
    int days_in_month[] = {31, 28 + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)), 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Calculate the month from the day of the year
    while (dayFraction > days_in_month[month-1]) {
        dayFraction -= days_in_month[month-1];
        month++;
    }

    // Calculate day fraction as month of day and fractional day from the part of day fraction < 1 and day_of_year (now day of month)
    //day_of_year = day_of_year + days_in_month[month-1]; // 1 month more gets removed above than necesserry
    day_of_year += dayFraction - floor(dayFraction);
    month++;

    // Math for Julian Date calculation
    if (month <= 2) {
        year = year - 1;
        month = month + 12;
    }
    B = 2 - floor( (double) year/100) + floor( (double) year / 400);
    JD = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + dayFraction + B - 1524.5;
    return JD;
}

double computeGMST(double jd) {
    double theta_g_0, theta_g, TU;
    double jd_0_UTC;

    // calculateion for 00:00 UTC
    jd_0_UTC = int(jd) + 0.5;
    TU = (jd_0_UTC - 2451545) / 36525;
    theta_g_0 = 24110.54841 + 8640184.812866 * TU + 0.093104 * TU * TU + 6.2 * pow(10,-6) * TU * TU * TU;

    // calculation for actual time
    double w_e, T, day_fraction;
    // to get the time in seconds since 00:00 we transform the day fractions to seconds
    day_fraction = jd - int(jd) - 0.5; 
    if (day_fraction < 0) day_fraction += 1;
    T = day_fraction * 24 * 60 * 60; 

    w_e = 7.292115 * pow (10, -5);
    theta_g = theta_g_0 + w_e * T;
    while (theta_g > TWO_PI) {
        theta_g -= TWO_PI;
    }
    while (theta_g <= 0) {
        theta_g += TWO_PI;
    }

    return theta_g;

}


