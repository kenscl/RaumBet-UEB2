#include "matlib.h"
#include <cstdio>
#include <cstdlib>

Vector_3D::Vector_3D() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
}

Vector_3D::Vector_3D(double x, double y, double z){
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector_3D::Vector_3D(Vector_3D const &other){
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
}

Vector_3D Vector_3D::cross_product(Vector_3D other) const{
    Vector_3D result;
    result.x = this->y * other.z - this->z * other.y;
    result.y = this->z * other.x - this->x * other.z;
    result.z = this->x * other.y - this->y * other.x;
    return result;
}

double Vector_3D::norm() const{
    double norm;
    norm = sqrt(x*x + y*y + z*z);
    return norm;
}
Vector_3D Vector_3D::normalize() const{
    double norm;
    norm = this->norm();
    if (norm == 0) return Vector_3D(0,0,0);

    Vector_3D result;
    result.x = this->x / norm;
    result.y = this->y / norm;
    result.z = this->z / norm;
    
    return result;
}

double Vector_3D::calculate_angle(Vector_3D other) const{
    double ac;
    Vector_3D v1, v2;
    v1 = this->normalize();
    v2 = other.normalize();
    double res;
    res = v1 * v2;
    ac = acos (res);
    return ac;
}

Matrix_3D Vector_3D::rmat_to(Vector_3D target) const{
    Vector_3D v1 = this->normalize();
    Vector_3D v2 = target.normalize();
    if (v1 == v2) {
        Matrix_3D ret;
        return ret.I();
    }
    Vector_3D vdiff = v2;

    double cos = v1 * v2;
    double length = cos;
    double sin = (v1.cross_product(v2)).norm();


    Matrix_3D rot;
    rot.r[0][0] = cos;
    rot.r[0][1] = - sin;
    rot.r[0][2] = 0;

    rot.r[1][0] = sin;
    rot.r[1][1] = cos;
    rot.r[1][2] = 0;

    rot.r[2][0] = 0; 
    rot.r[2][1] = 0; 
    rot.r[2][2] = 1; 

    vdiff.x -= length * v1.x;
    vdiff.y -= length * v1.y;
    vdiff.z -= length * v1.z;

    vdiff = vdiff.normalize();
    Vector_3D w = v2.cross_product(v1);

    Matrix_3D base_change_inv;

    base_change_inv.r[0][0] = v1.x;
    base_change_inv.r[0][1] = v1.y;
    base_change_inv.r[0][2] = v1.z;

    base_change_inv.r[1][0] = vdiff.x;
    base_change_inv.r[1][1] = vdiff.y;
    base_change_inv.r[1][2] = vdiff.z;

    base_change_inv.r[2][0] = w.x;
    base_change_inv.r[2][1] = w.y;
    base_change_inv.r[2][2] = w.z;

    Matrix_3D base_change;
    base_change = base_change_inv.inverse();
    Matrix_3D rmat = base_change * rot * base_change_inv;

    return rmat;

}

void Vector_3D::print() const{
    printf("Vector: \n");
    printf("x: %.5f \n", this->x);
    printf("y: %.5f \n", this->y);
    printf("z: %.5f \n", this->z);
}
Vector_3D Vector_3D::operator*(double value) const{
    Vector_3D result;
    result.x = value * this->x;
    result.y = value * this->y;
    result.z = value * this->z;
    return result;
}

Vector_3D Vector_3D::operator/(double value) const{
    Vector_3D result;
    result.x = this->x / value;
    result.y = this->y / value;
    result.z = this->z / value;
    return result;
}

double Vector_3D::operator*(Vector_3D other) const{
    double res;
    res = this->x * other.x + this->y * other.y + this->z * other.z;
    return res;
}

Vector_3D Vector_3D::operator+(Vector_3D other) const{
    Vector_3D result;
    result.x = this->x + other.x;
    result.y = this->y + other.y;
    result.z = this->z + other.z;
    return result;
}

Vector_3D Vector_3D::operator-(Vector_3D other) const{
    Vector_3D result;
    result.x = this->x - other.x;
    result.y = this->y - other.y;
    result.z = this->z - other.z;
    return result;
}

bool Vector_3D::operator==(Vector_3D other) const{
    if (this->x == other.x && this->y == other.y && this->z == other.z) return true;
    return false;
}

Matrix_3D::Matrix_3D(){
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            this->r[i][j] = 0;
        }
    }
}


Matrix_3D Matrix_3D::I() const{
    Matrix_3D ret;
    ret.r[0][0] = 1;
    ret.r[0][1] = 0;
    ret.r[0][2] = 0;

    ret.r[1][0] = 0;
    ret.r[1][1] = 1;
    ret.r[1][2] = 0;

    ret.r[2][0] = 0;
    ret.r[2][1] = 0;
    ret.r[2][2] = 1;
    return ret;
}
Vector_3D Matrix_3D::operator*(Vector_3D v) const{
    Vector_3D result;
    result.x = this->r[0][0] * v.x + this->r[0][1] * v.y + this->r[0][2] * v.z;
    result.y = this->r[1][0] * v.x + this->r[1][1] * v.y + this->r[1][2] * v.z;
    result.z = this->r[2][0] * v.x + this->r[2][1] * v.y + this->r[2][2] * v.z;
    return result;
}
Matrix_3D& Matrix_3D::operator=(const Matrix_3D& other){
        if (this != &other) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    r[i][j] = other.r[i][j];
                }
            }
        }
        return *this;
}
Matrix_3D Matrix_3D::operator*(const Matrix_3D& other) const{
    Matrix_3D result;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                result.r[i][j] += this->r[i][k] * other.r[k][j];
            }
        }
    }
    return result;
}

Matrix_3D Matrix_3D::operator*(const double scalar) const {
    Matrix_3D result;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result.r[i][j] = r[i][j] * scalar;
        }
    }
    return result;
}
Matrix_3D Matrix_3D::operator+(const Matrix_3D& other) const {
        Matrix_3D result;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result.r[i][j] = r[i][j] + other.r[i][j];
            }
        }
        return result;
    }
void Matrix_3D::print() const{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", r[i][j]);
            if (j == 2) printf("\n");
        }
    }
}

double Matrix_3D::determinant() const{
    double det = this->r[0][0] * (this->r[1][1] * this->r[2][2] - this->r[1][2] * this->r[2][1]) -
        this->r[0][1] * (this->r[1][0] * this->r[2][2] - this->r[1][2] * this->r[2][0]) +
        this->r[0][2] * (this->r[1][0] * this->r[2][1] - this->r[1][1] * this->r[2][0]);
    return det;
}
Matrix_3D Matrix_3D::inverse() const{
        double det = this->determinant();

        if (det == 0.0) {
            printf("matrix is singular no inversion possible \n");
        }

        Matrix_3D inv;

        inv.r[0][0] = (this->r[1][1] * this->r[2][2] - this->r[1][2] * this->r[2][1]) / det;
        inv.r[0][1] = (this->r[0][2] * this->r[2][1] - this->r[0][1] * this->r[2][2]) / det;
        inv.r[0][2] = (this->r[0][1] * this->r[1][2] - this->r[0][2] * this->r[1][1]) / det;

        inv.r[1][0] = (this->r[1][2] * this->r[2][0] - this->r[1][0] * this->r[2][2]) / det;
        inv.r[1][1] = (this->r[0][0] * this->r[2][2] - this->r[0][2] * this->r[2][0]) / det;
        inv.r[1][2] = (this->r[0][2] * this->r[1][0] - this->r[0][0] * this->r[1][2]) / det;

        inv.r[2][0] = (this->r[1][0] * this->r[2][1] - this->r[1][1] * this->r[2][0]) / det;
        inv.r[2][1] = (this->r[0][1] * this->r[2][0] - this->r[0][0] * this->r[2][1]) / det;
        inv.r[2][2] = (this->r[0][0] * this->r[1][1] - this->r[0][1] * this->r[1][0]) / det;

        return inv;
}

double lerp(double lly, double lry, double llx, double lrx, double x) {
    double distance = abs (lrx - llx);
    if (distance == 0) {
        return (lly + lry) / 2; 
    }
    double weight_left = lly * abs (lrx - x) / distance;
    double weight_right = lry * abs (llx - x) / distance;
    return weight_left + weight_right;
}

double linear_interpolation(double * x, double * y, double value, int len){
    if (x[0] > value) return y[0];
    for (int i = 0; i < len - 1; i++) {
        if (x[i] <= value && x[i+1] >= value) return lerp(y[i], y[i+1], x[i], x[i + 1], value);
    }
    return y[len];
}
