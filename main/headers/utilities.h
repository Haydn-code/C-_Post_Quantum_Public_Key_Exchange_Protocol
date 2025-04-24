#ifndef UTILITIES_H
#define UTILITIES_H

#include <array>
#include <stdexcept>
#include <iostream>
#include "montgomery_isogney.h"

class TorsionData {
	// class to generate torsion data
private:
	static void batch_cofactor_mul_generic(std::vector<Point>& G_list, fq_ctx_t& field, std::vector<slong>& pis, fmpz_t lower, fmpz_t upper);
	static void batch_cofactor_mul_generic(std::vector<fq_poly_t*>& G_list, fq_ctx_t& field, std::vector<slong>& pis, fmpz_t lower, fmpz_t upper);
	static std::vector<slong> has_order_constants(fmpz_t& D);
	static bool has_order_D(Point G, fmpz_t& D);
	static bool has_order_D(fq_poly_t& G, fq_ctx_t& field, fmpz_t& D);
	static Point generate_random_point(EllipticCurve& E, flint_rand_t& prng);
	static Point generate_point_order_D(EllipticCurve& E, fmpz_t& D, fmpz_t& p, flint_rand_t& prng);
	static Point generate_linearly_independant(EllipticCurve& E, Point& P, fmpz_t& D, fmpz_t& p, flint_rand_t& prng);
public:
	static std::array<Point, 2> torsion_basis(EllipticCurve& E, fmpz_t& D, fmpz_t& characteristic, flint_rand_t& prng);
};

#endif