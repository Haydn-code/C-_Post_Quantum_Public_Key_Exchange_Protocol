#ifndef MONTGOMERY_XZ_H
#define MONTGOMERY_XZ_H

#include <array>
#include <vector>

#include "elliptic_curve.h"

class KummerLine {
	// class for Kummer lines over montgomery curves
public:
	fq_ctx_t* field;
	fq_poly_t k_a;
	fq_poly_t k_c;
	KummerLine() { fmpz_t two; fmpz_set_ui(two, 2); fq_ctx_init(*field, two, 2, "x"); fq_poly_init(k_a, *field); fq_poly_init(k_c, *field); };
	KummerLine(fq_ctx_t& ring, std::array<int, 2> coeff);
	KummerLine(fq_ctx_t& ring, std::array<fq_poly_t*, 2> coeff);
	KummerLine(EllipticCurve& E);
	bool is_equal(KummerLine& other);
	fq_poly_t* j_invariant();
	fq_poly_t* a();
	fq_poly_t* curve_j_invariant();
	bool is_isomorphic(KummerLine& other);
};

class KummerPoint {
	// Class for Kummer points over montgomery curves
public:
	KummerLine p_parent;
	fq_poly_t p_x;
	fq_poly_t p_z;
	KummerPoint() { p_parent = KummerLine(); fq_poly_init(p_x, *p_parent.field); fq_poly_init(p_z, *p_parent.field); };
	KummerPoint(KummerLine& parent, fq_poly_t& x);
	KummerPoint(KummerLine& parent, std::array<fq_poly_t*, 2> coords);
	bool is_equal(KummerPoint& other);
	fq_poly_t* x();
	int curve_point();
	std::array<fq_poly_t*, 2> xDBL(fq_poly_t& X, fq_poly_t& Z, fq_poly_t& A, fq_poly_t& C);
	std::array<fq_poly_t*, 2> xADD(fq_poly_t XP, fq_poly_t ZP, fq_poly_t& XQ, fq_poly_t& ZQ, fq_poly_t& xPQ, fq_poly_t& zPQ);
	std::array<fq_poly_t*, 4> xDBLADD(fq_poly_t& XP, fq_poly_t& ZP, fq_poly_t& XQ, fq_poly_t& ZQ, fq_poly_t& xPQ, fq_poly_t& zPQ, fq_poly_t& A24, fq_poly_t& C24);
	KummerPoint twice();
	KummerPoint add(KummerPoint Q, KummerPoint PQ);
	KummerPoint multiply(fmpz_t m);
	KummerPoint ladder_3_pt(KummerPoint& xP, KummerPoint& xPQ, fmpz_t m);
	std::vector<KummerPoint> multiples();
};
#endif