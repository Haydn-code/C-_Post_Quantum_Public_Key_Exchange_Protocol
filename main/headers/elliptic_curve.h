#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include <array>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <iostream>
#include <flint/flint.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fmpz-conversions.h>

class EllipticCurve {
	// Class for an elliptic curve over a finite field
public:
	std::array<int, 5> curve_coefficients;
	fq_ctx_t* field;
	fmpz_t order;
	fmpz_t characteristic;
	EllipticCurve() { curve_coefficients = { 0, 6, 0, 1, 0 }; fmpz_set_ui(order, 1); fmpz_set_ui(characteristic, 3);  fq_ctx_init(*field, characteristic, 1, "i"); std::cout << "test" << std::endl; };
	EllipticCurve(fq_ctx_t& ring, std::array<int, 5> coefficients, fmpz_t& ord, fmpz_t& charact);
	bool equal(EllipticCurve other);
	~EllipticCurve();
};

class Point {
	// Class to represent a point on an elliptic curve
public:
	fq_poly_t p_x;
	fq_poly_t p_y;
	fq_poly_t p_z;
	fmpz_t order;
	EllipticCurve curve;
	Point() { std::cout << "called" << std::endl;  curve = EllipticCurve(); fq_poly_init(p_x, *curve.field); fq_poly_init(p_y, *curve.field); fq_poly_init(p_z, *curve.field); fmpz_init(order);}
	Point(fq_poly_t& i_x, fq_poly_t& i_y, fq_poly_t& i_z, EllipticCurve& p_curve);
	bool equal(Point other);
	static Point generate_random_element(EllipticCurve E, flint_rand_t& prng);
	static fq_poly_t* line_equation(EllipticCurve& E, Point& P, Point& Q, Point& D);
	static fq_poly_t* millers_loop(EllipticCurve& E, Point& P, Point Q, fmpz_t& D);
	static fq_poly_t* weil_pairing(EllipticCurve& E, Point& P, Point& Q, fmpz_t& D, flint_rand_t& prng);
	static fq_poly_t* line_slope(Point& P, Point& Q);
	static Point add_points(Point& P, Point& Q);
	static Point subtract(Point& P, Point& Q);
	static Point scalar_multiplication(Point& P, fmpz_t& m);
	~Point();
	fq_poly_t* y();
	fq_poly_t* x();
};

#endif