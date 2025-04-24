#include "../headers/elliptic_curve.h"

EllipticCurve::EllipticCurve(fq_ctx_t& ring, std::array<int, 5> coefficients, fmpz_t& ord, fmpz_t& charact) {
	// constructor for an elliptic curve

	curve_coefficients = coefficients;
	field = &ring;
	fmpz_set(order, ord);
	fmpz_set(characteristic, charact);
}

bool EllipticCurve::equal(EllipticCurve other) {
	// checks if an elliptic curve is equal to another elliptic curve

	if (std::equal(curve_coefficients.begin(), curve_coefficients.end(), other.curve_coefficients.begin())) {
		if (fmpz_equal(order, other.order)) {
			if (fmpz_equal(characteristic, other.characteristic)) {
				return true;
			}
		}
	}
	return false;
}

EllipticCurve::~EllipticCurve() {
	// destructor for an elliptic curve

	fmpz_clear(order);
	fmpz_clear(characteristic);
}

Point::Point(fq_poly_t& i_x, fq_poly_t& i_y, fq_poly_t& i_z, EllipticCurve& p_curve) {
	// constructor for a point on an elliptic curve (memory access errors on call)

	std::cout << "called" << std::endl;
	curve = p_curve;
	fq_poly_set(p_x, i_x, *curve.field);
	fq_poly_set(p_y, i_y, *curve.field);
	fq_poly_set(p_z, i_z, *curve.field);
	fmpz_set_ui(order, 1);
}

bool Point::equal(Point other) {
	// checks if a point is equal to another point

	if (curve.equal(other.curve)) {
		if (fq_poly_equal(p_x, other.p_x, *curve.field)) {
			if (fq_poly_equal(p_y, other.p_x, *curve.field)) {
				if (fq_poly_equal(p_z, other.p_x, *curve.field)) {
					return true;
				}
			}
		}
	}
	return false;
}

Point::~Point() {
	// destructor for a point on an elliptic curve

	fq_poly_clear(p_x, *curve.field);
	fq_poly_clear(p_y, *curve.field);
	fq_poly_clear(p_z, *curve.field);
	fmpz_clear(order);
	delete &curve;
}

fq_poly_t* Point::x() {
	// returns the affine x coordinate

	fq_poly_t q, r;
	fq_poly_divrem(q, r, p_x, p_z, *curve.field);
	fq_poly_add(q, q, r, *curve.field);
	fq_poly_clear(r, *curve.field);
	return &q;
}

fq_poly_t* Point::y() {
	// returns the affine y coordinate

	fq_poly_t q, r;
	fq_poly_divrem(q, r, p_y, p_z, *curve.field);
	fq_poly_add(q, q, r, *curve.field);
	fq_poly_clear(r, *curve.field);
	return &q;
}


Point Point::generate_random_element(EllipticCurve E, flint_rand_t& prng) {
	// generates a random element on curve with equation y^2 = x^3 + 6x^2 + x

	fq_poly_t x, x2, x3, y, z;
	fq_t multiplier;
	fq_poly_init(x, *E.field);
	fq_poly_init(x2, *E.field);
	fq_poly_init(x3, *E.field);
	fq_poly_init(y, *E.field);
	fq_poly_init(z, *E.field);
	fq_init(multiplier, *E.field);
	fq_poly_one(z, *E.field);

	std::cout << std::endl << "looping" << std::endl;
	// loops until value for x generated that has a corresponding y coordinate
	do {
		fq_poly_zero(y, *E.field); 
		fq_set_ui(multiplier, 6, *E.field);
		fq_poly_randtest(x, prng, 2, *E.field);
		fq_poly_pow(x2, x, 2, *E.field);
		fq_poly_scalar_mul_fq(x2, x2, multiplier, *E.field);
		fq_poly_pow(x3, x, 3, *E.field);
		fq_poly_add(y, x, x2, *E.field);
		fq_poly_add(y, y, x3, *E.field); 
	} while (!fq_poly_sqrt(y, y, *E.field));
	std::cout << "making point" << std::endl;
	Point p = Point(x, y, z, E);
	std::cout << std::endl << "made point" << std::endl;
	// deallocate memory
	fq_poly_clear(x, *E.field);
	fq_poly_clear(x2, *E.field);
	fq_poly_clear(x3, *E.field);
	fq_poly_clear(y, *E.field);
	fq_poly_clear(z, *E.field);
	std::cout << "returning" << std::endl;
	return p;
}

fq_poly_t* Point::line_equation(EllipticCurve& E, Point& P, Point& Q, Point& D) {
	// line equation algorithm for weil pairing sourced and adapted from https://github.com/Stetsyk/Weil-pairing-through-Miller-algorithm

	fq_poly_t result;
	if (fq_poly_equal(P.p_x, Q.p_x, *E.field)){
		if (fq_poly_equal(P.p_y, Q.p_y, *E.field)){
			if (fq_poly_is_zero(P.p_y, *E.field)) {
				fq_poly_sub(result, D.p_x, P.p_x, *E.field);
				return &result;
			}
		}
		else {
			fq_poly_sub(result, D.p_x, P.p_x, *E.field);
			return &result;
		}
	}
	fq_poly_t* lambda = line_slope(P, Q);
	fq_poly_t p, q;
	
	// numerator
	fq_poly_sub(p, D.p_x, P.p_x, *E.field);
	fq_poly_mul(p, *lambda, p, *E.field);
	fq_poly_sub(q, D.p_y, P.p_y, *E.field);
	fq_poly_sub(p, q, p, *E.field);

	// denominator
	fq_poly_mul(*lambda, *lambda, *lambda, *E.field);
	fq_poly_add(q, D.p_x, P.p_x, *E.field);
	fq_poly_add(q, q, Q.p_x, *E.field);
	fq_poly_sub(q, q, *lambda, *E.field);

	// result and memory cleanup
	fq_poly_divrem(result, q, p, q, *E.field);
	fq_poly_clear(*lambda, *E.field);
	fq_poly_clear(p, *E.field);
	fq_poly_clear(q, *E.field);

	return &result;
}

fq_poly_t* Point::millers_loop(EllipticCurve& E, Point& P, Point Q, fmpz_t& D) {
	// millers loop algorithm for weil pairing sourced and adapted from https://github.com/Stetsyk/Weil-pairing-through-Miller-algorithm

	fq_poly_t result;
	fq_poly_t* temp;
	fq_poly_one(result, *E.field);

	if (fq_poly_is_zero(P.p_z, *E.field) or fq_poly_is_zero(Q.p_z, *E.field) or (fq_poly_equal(P.p_y, Q.p_y, *E.field) and fq_poly_equal(P.p_x, Q.p_x, *E.field) and fq_poly_equal(P.p_z, Q.p_z, *E.field))) {
		return &result;
	}
	Point T = P;
	bool flag = false;
	// Iterate over the binary representation of the order
	char *bit_str = fmpz_get_str(NULL, 2, D);
	size_t i = strlen(bit_str);
	while (i >= 0) {
		// skips leading 0s
		if (!flag) {
			if (bit_str[i] == '1') {
				flag = true;
			}
			continue;
		}
		fq_poly_mul(result, result, result, *E.field);
		temp = line_equation(E, T, T, Q);
		fq_poly_mul(result, result, *temp, *E.field);
		T = add_points(T, T);
		if (bit_str[0] == 1) {
			temp = line_equation(E, T, P, Q);
			fq_poly_mul(result, result, *temp, *E.field);
			T = add_points(T, P);
		}
		--i;
	}

	return &result;
}

fq_poly_t* Point::weil_pairing(EllipticCurve& E, Point& P, Point& Q, fmpz_t& D, flint_rand_t& prng) {
	// weil pairing algorithm sourced and adapted from https://github.com/Stetsyk/Weil-pairing-through-Miller-algorithm

	Point S = generate_random_element(E, prng);
	fq_poly_t* fpqs;
	fq_poly_t* fps;
	fq_poly_t* fqps;
	fq_poly_t* fqs;

	// the functions for the points
	fpqs = millers_loop(E, P, add_points(Q, S), D);
	fps = millers_loop(E, P, S, D);
	fqps = millers_loop(E, P, subtract(P, S), D);
	fq_poly_neg(S.p_y, S.p_y, *S.curve.field);
	fqs = millers_loop(E, P, S, D);

	// calculate weil pairing
	fq_poly_mul(*fpqs, *fpqs, *fqs, *E.field);
	fq_poly_divrem(*fpqs, *fqs, *fpqs, *fps, *E.field);
	fq_poly_divrem(*fpqs, *fqs, *fpqs, *fqps, *E.field);

	// memory cleanup
	fq_poly_clear(*fps, *E.field);
	fq_poly_clear(*fqs, *E.field);
	fq_poly_clear(*fqps, *E.field);
	delete &S;

	return fpqs;
}

fq_poly_t* Point::line_slope(Point& P, Point& Q) {
	// calculates gradient of a line at a point or between two points

	fq_poly_t* P_x = P.x();
	fq_poly_t* P_y = P.y();
	fq_poly_t* Q_x = Q.x();
	fq_poly_t* Q_y = Q.y();
	fq_poly_t q;

	if (fq_poly_equal(*P_x, *Q_x, *P.curve.field) && fq_poly_equal(*P_y, *Q_y, *P.curve.field)) {
		// Point doubling (uses differentiation)
		fq_poly_t x, x2, y, r;
		fq_t multiplier;
		fq_poly_set(x, *P_x, *P.curve.field);
		fq_poly_pow(x2, x, 2, *P.curve.field);
		fq_set_ui(multiplier, 3, *P.curve.field);
		fq_poly_scalar_mul_fq(x2, x2, multiplier, *P.curve.field);
		fq_set_ui(multiplier, 12, *P.curve.field);
		fq_poly_scalar_mul_fq(x, x, multiplier, *P.curve.field);
		fq_set_ui(multiplier, 2, *P.curve.field);
		fq_poly_scalar_mul_fq(y, y, multiplier, *P.curve.field);
		fq_poly_add(x, x, x2, *P.curve.field);
		fq_poly_add_si(x, x, 1, *P.curve.field);
		fq_poly_divrem(q, r, x, y, *P.curve.field);
		fq_poly_clear(x, *P.curve.field);
		fq_poly_clear(x2, *P.curve.field);
		fq_poly_clear(y, *P.curve.field);
		fq_poly_clear(r, *P.curve.field);
	}
	else {
		// Point addition (calculates gradient of line between two points)
		fq_poly_t x_diff, y_diff, r;
		fq_poly_sub(x_diff, *P_x, *Q_x, *P.curve.field);
		fq_poly_sub(y_diff, *P_y, *Q_y, *P.curve.field);
		fq_poly_divrem(q, r, y_diff, x_diff, *P.curve.field);
		fq_poly_clear(x_diff, *P.curve.field);
		fq_poly_clear(y_diff, *P.curve.field);
		fq_poly_clear(r, *P.curve.field);
	}
	// deallocate memory
	fq_poly_clear(*P_x, *P.curve.field);
	fq_poly_clear(*P_y, *P.curve.field);
	fq_poly_clear(*Q_x, *P.curve.field);
	fq_poly_clear(*Q_y, *P.curve.field);

	return &q;
}

Point Point::add_points(Point& P, Point& Q) {
	// adds two point on an elliptic curve

	// any point plus the point at infinity is just itself
	if (fq_poly_is_zero(P.p_z, *P.curve.field)) {
		return Q;
	}

	if (fq_poly_is_zero(Q.p_z, *Q.curve.field)) {
		return P;
	}

	fq_poly_t* P_x = P.x();
	fq_poly_t* P_y = P.y();
	fq_poly_t* Q_x = Q.x();
	fq_poly_t* Q_y = Q.y();

	if (fq_poly_equal(*P_x, *Q_x, *P.curve.field)) {
		// checks if Q is -P
		fq_poly_t test;
		fq_poly_neg(test, *Q_x, *P.curve.field);
		if (fq_poly_equal(*P_x, test, *P.curve.field)) {
			// creates point at infinity
			fq_poly_clear(test, *P.curve.field);
			fq_poly_clear(*P_x, *P.curve.field);
			fq_poly_clear(*P_y, *P.curve.field);
			fq_poly_clear(*Q_x, *Q.curve.field);
			fq_poly_clear(*Q_y, *Q.curve.field);
			fq_poly_zero(P.p_x, *P.curve.field);
			fq_poly_one(P.p_y, *P.curve.field);
			fq_poly_zero(P.p_z, *P.curve.field);
			return P;
		}
	}

	fq_poly_t lambda2, x, y, z;
	fq_poly_t* lambda = line_slope(Q, P);

	// x coordinate
	fq_poly_pow(lambda2, *lambda, 2, *P.curve.field);
	fq_poly_sub(x, lambda2, *P_x, *P.curve.field);
	fq_poly_sub(x, x, *Q_x, *P.curve.field);

	// y coordinate
	fq_poly_sub(y, *P_x, x, *P.curve.field);
	fq_poly_mul(y, *lambda, y, *P.curve.field);
	fq_poly_sub(y, y, *P_y, *P.curve.field);

	// z coordinate
	fq_poly_one(z, *P.curve.field);

	// memory cleanup
	fq_poly_clear(lambda2, *P.curve.field);
	fq_poly_clear(*lambda, *P.curve.field);
	fq_poly_clear(*P_x, *P.curve.field);
	fq_poly_clear(*P_y, *P.curve.field);
	fq_poly_clear(*Q_x, *Q.curve.field);
	fq_poly_clear(*Q_y, *Q.curve.field);

	return Point(x, y, z, P.curve);
}

Point Point::subtract(Point& P, Point& Q) {
	// performs point subtraction on point P by Q

	Point inv_Q = Point(Q.p_x, Q.p_y, Q.p_z, Q.curve);
	fq_poly_neg(inv_Q.p_y, inv_Q.p_y, *inv_Q.curve.field);
	return add_points(P, inv_Q);
}

Point Point::scalar_multiplication(Point& P, fmpz_t& m) {
	// performs scalar multiplcation on point P by scalar m

	Point result = P;
	fmpz_t i;
	for (fmpz_zero(i); fmpz_cmp(i, m) > 0; fmpz_add_ui(i, i, 1)) {
		result = add_points(result, P);
	}
	fmpz_clear(i);
	fmpz_clear(m);
	return result;
}