#ifndef MONTGOMERY_ISOGENY_H
#define MONTGOMERY_ISOGENY_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "product_tree.h"
#include "montgomery_xz.h"

class Helper {
	// Class to provide helper functions throughout the code
public:
	static size_t find_index(std::vector<fmpz_t*> factor, fmpz_t search);
	static std::vector<std::array<fmpz_t*, 2>> factoring(fmpz_t degree);
};

class KummerLineIsogenyGeneric {
	// Parent class for all variants of Kummer line isogenies
public:
	fmpz_t i_degree;
	KummerPoint i_kernel;
	KummerLine i_domain;
	KummerLine i_codomain;
	KummerLineIsogenyGeneric() { fmpz_init(i_degree); i_kernel = KummerPoint(); i_domain = KummerLine(); i_codomain = KummerLine(); };
	virtual KummerPoint call(KummerPoint& P) { return KummerPoint(); };
};

class EvaluateIsogenies {
	// Class to evaluate isogenies
public:
	static KummerPoint evaluate_factored_kummer_isogeny(std::vector<KummerLineIsogenyGeneric*> phi_list, KummerPoint P);
	static std::vector<KummerLineIsogenyGeneric*> factored_kummer_isogeny(KummerLine K, KummerPoint P, fmpz_t& order, int threshold = 1000);
	static std::vector<KummerLineIsogenyGeneric*> sparse_isogeny_prime_power(KummerPoint P, fmpz_t& l, fmpz_t& e, float split = 0.8, int threshold = 1000);
	static std::vector<KummerLineIsogenyGeneric*> recursive_sparse_isogeny(KummerPoint& Q, fmpz_t l, fmpz_t k, float split, bool& algorithm);
};

class KummerLineIsogenyComposite : public KummerLineIsogenyGeneric{
	// Class for storage of composite isogenies of Kummer lines
public:
	std::vector<KummerLineIsogenyGeneric*> i_phis;
	fmpz_t i_degree;
	KummerLine i_domain;
	KummerLine i_codomain;
	KummerLineIsogenyComposite(KummerLine& domain, KummerPoint& kernel, fmpz_t& degree, int threshold=1500);
	KummerPoint call(KummerPoint& P) override;

};

class KummerLineIsogeny : public KummerLineIsogenyGeneric{
	// Class for storage of regular isogenies of Kummer lines
public:
	fmpz_t i_degree;
	KummerPoint i_kernel;
	KummerLine i_domain;
	KummerLine i_codomain;
	std::vector<std::array<fq_poly_t*, 2>> i_edwards_multiples;
	KummerLineIsogeny(KummerLine& domain, KummerPoint& kernel, fmpz_t degree);
	KummerPoint call(KummerPoint& P) override;
	std::vector<std::array<fq_poly_t*, 2>> precompute_edwards_multiples(fmpz_t& d);
	std::array<fq_poly_t*, 2> compute_codomain_constants();
	std::array<fq_poly_t*, 2> compute_codomain_constants_even();
	KummerLine compute_codomain();
	KummerPoint evaluate_isogeny(KummerPoint& P);
	KummerPoint evaluate_isogeny_even(KummerPoint& P);
};

class KummerLineIsogeny_VeluSqrt : public KummerLineIsogenyGeneric{
	// Class for Kummer line isogenies using velu's formulas for speedup
public:
	fmpz_t i_degree;
	KummerPoint i_kernel;
	KummerLine i_domain;
	KummerLine i_codomain;
	fq_poly_t i_a;
	fq_poly_t i_ring_generator;
	fq_ctx_t* i_ring;
	fmpz_t stop;
	ProductTree hI_tree;
	std::vector<std::array<fq_poly_t*, 3>> EJ_parts;
	fq_poly_t* hK;
	KummerLineIsogeny_VeluSqrt(KummerLine& domain, KummerPoint& kernel, fmpz_t degree);
	KummerPoint call(KummerPoint& P) override { return evaluate_isogeny(P); };
	fmpz_t* hI_resultant(fq_poly_t& poly);
	ProductTree hI_precomputation(fmpz_t& b, fmpz_t& c);
	std::array<fq_poly_t*, 3> Fs(fq_poly_t& X1, fq_poly_t& X2);
	std::vector<std::array<fq_poly_t*, 3>> EJ_precomputation(fmpz_t& b);
	fq_poly_t* hK_precomputation(fmpz_t& b,fmpz_t& c);
	std::array<fq_poly_t*, 2> compute_codomain_constants();
	KummerLine compute_codomain();
	KummerPoint evaluate_isogeny(KummerPoint& P);
};
#endif