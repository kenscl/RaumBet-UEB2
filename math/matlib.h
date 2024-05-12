#ifndef __MATLIB
#define __MATLIB
# include "math.h"
# define PI                 3.14159265359
# define C                  299792458.0
# define PLANCK_CONSTANT    6.62607015e-34
# define BOLTZMANN_CONSTANT 1.380649e-23

class Vector_3D;
class Matrix_3D;
class Vector_3D {
    public: 
        double x,y,z;

    public: 
        Vector_3D();
        Vector_3D(double x, double y, double z);
        Vector_3D(Vector_3D const &other);

        Vector_3D cross_product(Vector_3D other) const;
        double norm() const;
        Vector_3D normalize() const;
        void print() const;
        double calculate_angle(Vector_3D other) const;
        Matrix_3D rmat_to(Vector_3D target) const;

        Vector_3D operator*(double value) const;
        Vector_3D operator/(double value) const;
        bool operator==(Vector_3D other) const;

        double operator*(Vector_3D other) const;
        Vector_3D operator+(Vector_3D other) const;
        Vector_3D operator-(Vector_3D other) const;
};
class Matrix_3D {
    public: 
        double r[3][3];

    public: 
        Matrix_3D();

        Matrix_3D inverse() const;
        double determinant() const;
        Matrix_3D I() const;
        void print() const;

        Vector_3D operator*(Vector_3D v) const;
        Matrix_3D& operator=(const Matrix_3D& other);
        Matrix_3D operator*(const Matrix_3D& other) const;
        Matrix_3D operator*(const double scalar) const;
        Matrix_3D operator+(const Matrix_3D& other) const;
        
};

double lerp(double lly, double lry, double llx, double lrx, double x);
double linear_interpolation(double * x, double * y, double value, int len);
#endif

