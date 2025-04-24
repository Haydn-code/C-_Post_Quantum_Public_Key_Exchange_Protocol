#include <bitset>
#include <vector>

#include "../headers/montgomery_xz.h"

KummerLine::KummerLine(fq_ctx_t& ring, std::array<int, 2> coeff) {
	// constructor for kummer line given a and c invariants as integers

	field = &ring;
	fq_t temp1, temp2;
	fq_set_ui(temp1, coeff[0], *field);
	fq_set_ui(temp2, coeff[1], *field);
	fq_poly_set_coeff(k_a, 0, temp1, *field);
	fq_poly_set_coeff(k_c, 0, temp2, *field);

	// deallocate memory
	fq_clear(temp1, *field);
	fq_clear(temp2, *field);
};

KummerLine::KummerLine(fq_ctx_t& ring, std::array<fq_poly_t*, 2> coeff) {
	// constructor for a kumerline given a and c invariants as polynomials

	field = &ring;
	fq_poly_set(k_a, *coeff[0], *field);
	fq_poly_set(k_c, *coeff[1], *field);
}

KummerLine::KummerLine(EllipticCurve& E) {
	// constructor for kummer line given an elliptic curve

	field = E.field; 
	fq_t temp1, temp2;
	fq_set_ui(temp1, E.curve_coefficients[1], *field);
	fq_set_ui(temp2, 1, *field);
	fq_poly_set_coeff(k_a, 1, temp1, *field);
	fq_poly_set_coeff(k_c, 1, temp2, *field);

	//deallocate dynamic memory
	fq_clear(temp1, *field);
	fq_clear(temp2, *field);
}

bool KummerLine::is_equal(KummerLine& other) {
	// checks whether a kummer line is equal to another kummer line

	fq_poly_t test1, test2;
	fq_poly_mul(test1, k_a, other.k_c, *field);
	fq_poly_mul(test2, k_c, other.k_a, *field);
	if (field == other.field && fq_poly_equal(test1, test2, *field)) {
		fq_poly_clear(test1, *field);
		fq_poly_clear(test2, *field);
		return true;
	}
	fq_poly_clear(test1, *field);
	fq_poly_clear(test2, *field);
	return false;
}

fq_poly_t* KummerLine::j_invariant() {
	// calculates the j invariant of a kummerline

	fq_t multiplier;
	fq_poly_t temp1, temp2, j_num, j_den, j;
	
	// calculate numerator
	fq_poly_pow(temp1, k_a, 2, *field);
	fq_poly_pow(temp2, k_c, 2, *field);
	fq_set_ui(multiplier, 3, *field);
	fq_poly_scalar_mul_fq(temp2, temp2, multiplier, *field);
	fq_poly_sub(temp1, temp1, temp2, *field);
	fq_poly_pow(temp1, temp1, 3, *field);
	fq_set_ui(multiplier, 256, *field);
	fq_poly_scalar_mul_fq(j_num, temp1, multiplier, *field);

	// calculate denomenator
	fq_poly_pow(temp1, k_a, 2, *field);
	fq_poly_pow(temp2, k_c, 2, *field);
	fq_poly_sub(temp1, temp1, temp2, *field);
	fq_poly_pow(temp2, k_c, 4, *field);
	fq_poly_mul(j_den, temp2, temp1, *field);

	// calculate j invariant
	fq_poly_divrem(j, temp1, j_num, j_den, *field);

	// deallocate dynamic memory
	fq_clear(multiplier, *field);
	fq_poly_clear(temp1, *field);
	fq_poly_clear(temp2, *field);
	fq_poly_clear(j_num, *field);
	fq_poly_clear(j_den, *field);

	return &j;
}

fq_poly_t* KummerLine::a() {
	// calculate a invariant

	fq_poly_t q, r;
	fq_poly_divrem(q, r, k_a, k_c, *field);
	fq_poly_clear(r, *field);
	return &q;
}

fq_poly_t* KummerLine::curve_j_invariant() {
	// calculates the j invariant of a Mongtomery elliptic curve according to sage maths implementation https://github.com/sagemath/sage/blob/develop/src/sage/schemes/elliptic_curves/ell_generic.py#L1146

	// coefficients of the elliptic curve created from a KummerLine
	fq_poly_t* coeff = a();
	fq_poly_t one;
	fq_poly_one(one, *field);

	fq_poly_t b2, b4, b8;
	fq_t multiplier;

	// calculate b invariants
	fq_set_ui(multiplier, 4, *field);
	fq_poly_scalar_mul_fq(b2, *coeff, multiplier, *field);

	fq_set_ui(multiplier, 2, *field);
	fq_poly_scalar_mul_fq(b4, one, multiplier, *field);

	fq_poly_pow(b8, one, 2, *field);
	fq_poly_neg(b8, b8, *field);

	// calculate denominator
	fq_poly_t den, temp;

	fq_poly_pow(den, b2, 2, *field);
	fq_poly_mul(den, den, b8, *field);
	fq_poly_neg(den, den, *field);

	fq_poly_pow(temp, b4, 3, *field);
	fq_set_ui(multiplier, 8, *field);
	fq_poly_scalar_mul_fq(temp, temp, multiplier, *field);
	fq_poly_sub(den, den, temp, *field);

	// calculate nominator
	fq_poly_t nom;

	fq_poly_pow(nom, b2, 2, *field);
	fq_set_ui(multiplier, 24, *field);
	fq_poly_scalar_mul_fq(temp, b4, multiplier, *field);
	fq_poly_sub(nom, nom, temp, *field);
	fq_poly_pow(nom, nom, 3, *field);

	// calculate j invariant
	fq_poly_divrem(nom, den, nom, den, *field);

	// deallocate unused dynamic memory
	fq_poly_clear(*coeff, *field);
	fq_poly_clear(one, *field);
	fq_poly_clear(b2, *field);
	fq_poly_clear(b4, *field);
	fq_poly_clear(b8, *field);
	fq_poly_clear(den, *field);
	fq_poly_clear(temp, *field);
	fq_clear(multiplier, *field);

	return &nom;
}

bool KummerLine::is_isomorphic(KummerLine& other) {
	// checks if kummerline is isomorphic with another kummerline

	return fq_poly_equal(*j_invariant(), *other.j_invariant(),*field);
}

KummerPoint::KummerPoint(KummerLine& parent, std::array<fq_poly_t*, 2> coords) {
	// constructor for a kummerpoint given a kummerline and two coordinates

	p_parent = parent;
	fq_poly_set(p_x, *coords[0], *p_parent.field);
	fq_poly_set(p_z, *coords[1], *p_parent.field);
}

KummerPoint::KummerPoint(KummerLine& parent, fq_poly_t& x) {
	// constructor for a kummerpoint given a kummerline and an affine x coordinate

	p_parent = parent;
	fq_poly_set(p_x, x, *p_parent.field);
	fq_poly_one(p_z, *p_parent.field);
}

bool KummerPoint::is_equal(KummerPoint& other) {
	// checks if a kummerpoint is equal to another kummerpoint

	fq_poly_t test1, test2;
	fq_poly_mul(test1, p_x, other.p_z, *p_parent.field);
	fq_poly_mul(test2, p_z, other.p_x, *p_parent.field);
	if (p_parent.field == other.p_parent.field && fq_poly_equal(test1, test2, *p_parent.field)) {
		fq_poly_clear(test1, *p_parent.field);
		fq_poly_clear(test2, *p_parent.field);
		return true;
	}
	fq_poly_clear(test1, *p_parent.field);
	fq_poly_clear(test2, *p_parent.field);
	return false;
}

fq_poly_t* KummerPoint::x() {
	// finds the affine x coordinate

	fq_poly_t q, r;
	fq_poly_divrem(q, r, p_x, p_z, *p_parent.field);
	fq_poly_clear(r, *p_parent.field);
	return &q;
}

std::array<fq_poly_t*, 2> KummerPoint::xDBL(fq_poly_t& X, fq_poly_t& Z, fq_poly_t& A, fq_poly_t& C) {
	// function for montgomery doubling of a point with projecitve coordinates (X:Z) with projective curve constants (A:C)

	fq_poly_t t0, t1, Z2, X2;
	// point doubling algorithm for x only montgomery
	fq_poly_sub(t0, X, Z, *p_parent.field);
	fq_poly_add(t0, X, Z, *p_parent.field);
	fq_poly_pow(t0, t0, 2, *p_parent.field);
	fq_poly_pow(t1, t1, 2, *p_parent.field);
	fq_poly_mul(Z2, C, t0, *p_parent.field);
	fq_poly_add(Z2, Z2, Z2, *p_parent.field);
	fq_poly_add(Z2, Z2, Z2, *p_parent.field);
	fq_poly_mul(X2, Z2, t1, *p_parent.field);
	fq_poly_sub(t1, t1, t0, *p_parent.field);
	fq_poly_add(t0, C, C, *p_parent.field);
	fq_poly_add(t0, t0, A, *p_parent.field);
	fq_poly_mul(t0, t0, t1, *p_parent.field);
	fq_poly_add(Z2, Z2, t0, *p_parent.field);
	fq_poly_mul(Z2, Z2, t1, *p_parent.field);

	// deallocate unused dynamic memory
	fq_poly_clear(t0, *p_parent.field);
	fq_poly_clear(t1, *p_parent.field);

	return { &X2, &Z2 };
}

std::array<fq_poly_t*, 2> KummerPoint::xADD(fq_poly_t XP, fq_poly_t ZP, fq_poly_t& XQ, fq_poly_t& ZQ, fq_poly_t& xPQ, fq_poly_t& zPQ) {
	// function for montgomery differential addition

	fq_poly_t t0, t1, XQP, ZQP;
	fq_poly_add(t0, XP, ZP, *p_parent.field);
	fq_poly_sub(t1, XP, ZP, *p_parent.field);
	fq_poly_sub(XP, XQ, ZQ, *p_parent.field);
	fq_poly_add(ZP, XQ, ZQ, *p_parent.field);
	fq_poly_mul(t0, XP, t0, *p_parent.field);
	fq_poly_mul(t1, ZP, t1, *p_parent.field);
	fq_poly_sub(ZP, t0, t1, *p_parent.field);
	fq_poly_add(XP, t0, t1, *p_parent.field);
	fq_poly_pow(ZP, ZP, 2, *p_parent.field);
	fq_poly_pow(XQP, XP, 2, *p_parent.field);
	fq_poly_mul(ZQP, xPQ, ZP, *p_parent.field);
	fq_poly_mul(XQP, XQP, zPQ, *p_parent.field);
	
	// deallocate unused dynamic memory
	fq_poly_clear(XP, *p_parent.field);
	fq_poly_clear(ZP, *p_parent.field);
	fq_poly_clear(t0, *p_parent.field);
	fq_poly_clear(t1, *p_parent.field);

	return { &XQP, &ZQP };
}

std::array<fq_poly_t*, 4> KummerPoint::xDBLADD(fq_poly_t& XP, fq_poly_t& ZP, fq_poly_t& XQ, fq_poly_t& ZQ, fq_poly_t& xPQ, fq_poly_t& zPQ, fq_poly_t& A24, fq_poly_t& C24) {
	// function for step in montgomery ladder, simultaneous doubling and differential addition

	fq_poly_t t0, t1, t2, X2P, XQP, Z2P, ZQP;
	fq_poly_add(t0, XP, ZP, *p_parent.field);
	fq_poly_sub(t1, XP, ZP, *p_parent.field);
	fq_poly_pow(X2P, t0, 2, *p_parent.field);
	fq_poly_sub(t2, XQ, ZQ, *p_parent.field);
	fq_poly_add(XQP, XQ, ZQ, *p_parent.field);
	fq_poly_mul(t0, t0, t2, *p_parent.field);
	fq_poly_pow(Z2P, t1, 2, *p_parent.field);
	fq_poly_mul(t1, t1, XQP, *p_parent.field);
	fq_poly_sub(t2, X2P, Z2P, *p_parent.field);
	fq_poly_mul(Z2P, Z2P, C24, *p_parent.field);
	fq_poly_mul(X2P, X2P, Z2P, *p_parent.field);
	fq_poly_mul(XQP, A24, t2, *p_parent.field);
	fq_poly_sub(ZQP, t0, t1, *p_parent.field);
	fq_poly_add(Z2P, XQP, Z2P, *p_parent.field);
	fq_poly_add(XQP, t0, t1, *p_parent.field);
	fq_poly_mul(Z2P, Z2P, t2, *p_parent.field);
	fq_poly_pow(ZQP, ZQP, 2, *p_parent.field);
	fq_poly_pow(XQP, XQP, 2, *p_parent.field);
	fq_poly_mul(ZQP, xPQ, ZQP, *p_parent.field);
	fq_poly_mul(XQP, XQP, zPQ, *p_parent.field);

	// deallocate unnecessary memory
	fq_poly_clear(t0, *p_parent.field);
	fq_poly_clear(t1, *p_parent.field);
	fq_poly_clear(t2, *p_parent.field);

	return { &X2P, &Z2P, &XQP, &ZQP };
}

KummerPoint KummerPoint::twice() {
	// doubles a kummerpoint

	if (not p_z) {
		return KummerPoint(p_parent, { &p_x, &p_z });
	}
	std::array<fq_poly_t*, 2> xz = xDBL(p_x, p_z, p_parent.k_a, p_parent.k_c);
	return KummerPoint(p_parent, xz);
}

KummerPoint KummerPoint::add(KummerPoint Q, KummerPoint PQ) {
	// adds two kummerpoints together

	if (not p_z) {
		return Q;
	}
	if (not Q.p_z) {
		return KummerPoint(p_parent, { &p_x, &p_z });
	}
	if (not PQ.p_z) {
		return twice();
	}
	std::array<fq_poly_t*, 2> xz = xADD(p_x, p_z, Q.p_x, Q.p_z, PQ.p_x, PQ.p_z);
	return KummerPoint(p_parent, xz);
}

KummerPoint KummerPoint::multiply(fmpz_t m) {
	// multiplies a kummerpoint by a scalar integer m

	if (fmpz_is_zero(m)) {
		throw std::invalid_argument("m cannot be equal to 0");
	}
	fq_poly_t XP, ZP, X0, X1, Z0, Z1, A, C, A24, C24;
	fmpz_abs(m, m);
	
	// initiate variables for loop
	fq_poly_set(XP, p_x, *p_parent.field);
	fq_poly_set(ZP, p_z, *p_parent.field);
	fq_poly_one(X0, *p_parent.field);
	fq_poly_one(Z0, *p_parent.field);
	fq_poly_set(X1, XP, *p_parent.field);
	fq_poly_set(Z1, ZP, *p_parent.field);

	// converting parameters for projective DBLADD -> (A24:C24)=(A+2C:4C)
	fq_poly_set(A, p_parent.k_a, *p_parent.field);
	fq_poly_set(C, p_parent.k_c, *p_parent.field);
	fq_poly_add(A24, C, C, *p_parent.field);
	fq_poly_add(C24, A24, A24, *p_parent.field);
	fq_poly_add(A24, A24, A, *p_parent.field);


	// montgomery ladder
	char* bit_str = fmpz_get_str(NULL, 2, m);
	for (size_t i = 0; i < strlen(bit_str); ++i) {
		if (bit_str[i] == '0') {
			std::array<fq_poly_t*, 4> temp = xDBLADD(X0, Z0, X1, Z1, XP, ZP, A24, C24);
			fq_poly_set(X0, *temp[0], *p_parent.field);
			fq_poly_set(Z0, *temp[1], *p_parent.field);
			fq_poly_set(X1, *temp[2], *p_parent.field);
			fq_poly_set(Z1, *temp[3], *p_parent.field);
		}
		else {
			std::array<fq_poly_t*, 4> temp = xDBLADD(X1, Z1, X0, Z0, XP, ZP, A24, C24);
			fq_poly_set(X1, *temp[0], *p_parent.field);
			fq_poly_set(Z1, *temp[1], *p_parent.field);
			fq_poly_set(X0, *temp[2], *p_parent.field);
			fq_poly_set(Z0, *temp[3], *p_parent.field);
		}
	}

	// deallocate unused dynamic memory
	fmpz_clear(m);
	fq_poly_clear(XP, *p_parent.field);
	fq_poly_clear(ZP, *p_parent.field);
	fq_poly_clear(X1, *p_parent.field);
	fq_poly_clear(Z1, *p_parent.field);
	fq_poly_clear(A, *p_parent.field);
	fq_poly_clear(C, *p_parent.field);
	fq_poly_clear(A24, *p_parent.field);
	fq_poly_clear(C24, *p_parent.field);

	return KummerPoint(p_parent, { &X0, &Z0 });
}

KummerPoint KummerPoint::ladder_3_pt(KummerPoint& xP, KummerPoint& xPQ, fmpz_t m) {
	// computes xP + [m]*xQ where self is xQ with x only arithmetic 

	if (fmpz_is_zero(m)) {
		return xP;
	}
	fmpz_abs(m, m);
	fq_poly_t A, C, A24, C24, XQ, ZQ, XP, ZP, XPQ, ZPQ;

	// converting parmeters for projective DBLADD -> (A24:C24)=(A+2C:4C)
	fq_poly_set(A, p_parent.k_a, *p_parent.field);
	fq_poly_set(C, p_parent.k_c, *p_parent.field);
	fq_poly_add(A24, C, C, *p_parent.field);
	fq_poly_add(C24, A24, A24, *p_parent.field);
	fq_poly_add(A24, A24, A, *p_parent.field);

	// initialise variables for 3 point ladder
	fq_poly_set(XQ, p_x, *p_parent.field);
	fq_poly_set(ZQ, p_z, *p_parent.field);
	fq_poly_set(XP, xP.p_x, *xP.p_parent.field);
	fq_poly_set(ZP, xP.p_z, *xP.p_parent.field);
	fq_poly_set(XPQ, xPQ.p_x, *xPQ.p_parent.field);
	fq_poly_set(ZPQ, xPQ.p_z, *xPQ.p_parent.field);
	
	// compute 3 point ladder
	char* bit_str = fmpz_get_str(NULL, 2, m);
	size_t i = strlen(bit_str);
	while (i >= 0) {
		if (bit_str[i] == '1') {
			std::array<fq_poly_t*, 4> temp = xDBLADD(XQ, ZQ, XP, ZP, XPQ, ZPQ, A24, C24);
			fq_poly_set(XQ, *temp[0], *p_parent.field);
			fq_poly_set(ZQ, *temp[1], *p_parent.field);
			fq_poly_set(XP, *temp[2], *xP.p_parent.field);
			fq_poly_set(ZP, *temp[3], *xP.p_parent.field);
		}
		else {
			std::array<fq_poly_t*, 4> temp = xDBLADD(XQ, ZQ, XPQ, ZPQ, XP, ZP, A24, C24);
			fq_poly_set(XQ, *temp[0], *p_parent.field);
			fq_poly_set(ZQ, *temp[1], *p_parent.field);
			fq_poly_set(XPQ, *temp[2], *xPQ.p_parent.field);
			fq_poly_set(ZPQ, *temp[3], *xPQ.p_parent.field);
		}
		--i;
	}

	// deallocate unused dynamic memory
	fmpz_clear(m);
	fq_poly_clear(A, *p_parent.field);
	fq_poly_clear(C, *p_parent.field);
	fq_poly_clear(A24, *p_parent.field);
	fq_poly_clear(C24, *p_parent.field);
	fq_poly_clear(XQ, *p_parent.field);
	fq_poly_clear(ZQ, *p_parent.field);
	fq_poly_clear(XPQ, *xPQ.p_parent.field);
	fq_poly_clear(ZPQ, *xPQ.p_parent.field);

	return KummerPoint(p_parent, { &XP, &ZP });
}

std::vector<KummerPoint> KummerPoint::multiples() {
	// returns all multiples of a KummerPoint

	std::vector<KummerPoint> multiples;
	KummerPoint self = KummerPoint(p_parent, { &p_x, &p_z });
	multiples.push_back(self);
	KummerPoint R_0 = twice();
	if (self.is_equal(R_0)) {
		return multiples;
	}
	KummerPoint Q = KummerPoint(p_parent, { &p_x, &p_z });
	while (true) {
		multiples.push_back(R_0);
		KummerPoint R_1 = R_0.add(self, Q);
		if (R_0.is_equal(R_1)) {
			return multiples;
		}
		Q = R_0;
		R_0 = R_1;
	}
}