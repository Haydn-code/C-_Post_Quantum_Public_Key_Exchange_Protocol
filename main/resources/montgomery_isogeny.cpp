#include "../headers/montgomery_isogney.h"

size_t Helper::find_index(std::vector<fmpz_t*> factor, fmpz_t search){
	// finds index of an element equal to an element associated at an address of a vector of addresses

	std::vector<fmpz_t*>::iterator start = factor.begin();
	std::vector<fmpz_t*>::iterator end = factor.end();
	size_t i = 0;
	while (start != end and fmpz_cmp(search, **start)) {
		++start;
		++i;
	}
	fmpz_clear(search);
	return size_t(i);
}

std::vector<std::array<fmpz_t*, 2>> Helper::factoring(fmpz_t degree) {
	// makes use of the trail prime factorisation algorithm to find the prime factor of degree

	std::vector<std::array<fmpz_t*, 2>> factors;
	if (fmpz_cmp_ui(degree, 1) <= 0) {
		return factors;
	}
	fmpz_t z;
	fmpz_t z2;
	fmpz_set_ui(z, 2);
	std::vector<fmpz_t*> factor;
	std::array<fmpz_t*, 2> temp;
	fmpz_set_ui(*temp[1], 1);

	// removes all 2 factors from degree
	while (fmpz_divisible(degree, z)) {
		size_t idx = find_index(factor, z);
		if (fmpz_cmp(*factor[idx], z) == 0) {
			fmpz_add_ui(*factors[idx][1], *factors[idx][1], 1);
		}
		else {
			fmpz_set(*temp[0], z);
			factors.push_back(temp);
			factor.push_back(&z);

		}
		fmpz_divexact(degree, degree, z);
	}

	fmpz_set_ui(z, 3);

	// trail prime factorisation algorithm
	fmpz_mul(z2, z, z);
	while (fmpz_cmp(z2, degree) <= 0){
		fmpz_mul(z2, z, z);
		if (fmpz_divisible(degree, z)) {
			size_t idx = find_index(factor, z);
			if (fmpz_cmp(*factor[idx], z) == 0) {
				fmpz_add_ui(*factors[idx][1], *factors[idx][1], 1);
			}
			else {
				fmpz_set(*temp[0], z);
				factors.push_back(temp);
				factor.push_back(&z);

			}
			fmpz_divexact(degree, degree, z);
		}
		else {
			fmpz_add_ui(z, z, 2);
			fmpz_mul(z2, z, z);
		}
	}
	
	// catches the last prime in the algorithm
	if (fmpz_cmp_ui(degree, 1) > 0) {
		size_t idx = find_index(factor, degree);
		if (fmpz_cmp(*factor[idx], degree) == 0) {
			fmpz_add_ui(*factors[idx][1], *factors[idx][1], 1);
		}
		else {
			fmpz_set(*temp[0], degree);
			factors.push_back(temp);
		}
	}

	// deallocates dynamic memory
	fmpz_clear(degree);
	fmpz_clear(z);
	fmpz_clear(z2);
	factor.clear();
	return factors;
}

KummerPoint EvaluateIsogenies::evaluate_factored_kummer_isogeny(std::vector<KummerLineIsogenyGeneric*> phi_list, KummerPoint P) {
	// evaluates a KummerPoint computed through a chain of isogenies

	for (int i = 0; i < phi_list.size(); i++) {
		P = phi_list[i]->call(P);
	}
	return P;
}

std::vector<KummerLineIsogenyGeneric*> EvaluateIsogenies::recursive_sparse_isogeny(KummerPoint& Q, fmpz_t l, fmpz_t k, float split, bool& algorithm) {
	// helper function for factored kummer isogeny 

	if (fmpz_is_one(k)) {
		if (algorithm) {
			KummerLineIsogeny_VeluSqrt* phi = new KummerLineIsogeny_VeluSqrt(Q.p_parent, Q, l);
			return { phi };
		}
		else {
			KummerLineIsogeny* phi = new KummerLineIsogeny(Q.p_parent, Q, l);
			return { phi };
		}
	}
	fmpz_t k1;
	fmpz_set_d(k1, (fmpz_get_d(k) * split + 0.5));
	fmpz_sub_ui(k, k, 1);
	if (fmpz_cmp(k, k1) < 0 or fmpz_cmp_ui(k1, 1) < 0) {
		if (fmpz_cmp_ui(k, 1) > 0) {
			fmpz_set(k1, k);
		}
		else {
			fmpz_one(k1);
		}
	}
	fmpz_add_ui(k, k, 1);
	fmpz_pow_fmpz(l, l, k1);
	KummerPoint Q1 = Q.multiply(l);
	fmpz_t param;
	fmpz_sub(param, k, k1);
	std::vector<KummerLineIsogenyGeneric*> L = recursive_sparse_isogeny(Q1, l, param, split, algorithm);

	KummerPoint Q2 = evaluate_factored_kummer_isogeny(L, Q);
	std::vector<KummerLineIsogenyGeneric*> R = recursive_sparse_isogeny(Q2, l, k1, split, algorithm);

	L.insert(L.end(), R.begin(), R.end());

	fmpz_clear(l);
	fmpz_clear(k);
	fmpz_clear(param);
	fmpz_clear(k1);

	return L;
}


std::vector<KummerLineIsogenyGeneric*> EvaluateIsogenies::sparse_isogeny_prime_power(KummerPoint P, fmpz_t& l, fmpz_t& e, float split, int threshold) {
	// determines the type of isogeny to initialise

	bool algorithm = false;
	if (fmpz_cmp_ui(l, threshold) > 0) {
		algorithm = true;
	}
	return recursive_sparse_isogeny(P, l, e, split, algorithm);
}

std::vector<KummerLineIsogenyGeneric*> EvaluateIsogenies::factored_kummer_isogeny(KummerLine K, KummerPoint P, fmpz_t& order, int threshold) {
	// returns a chain of isogenies generated from  the factors of a large order

	fmpz_t cofactor, temp;
	fmpz_set(cofactor, order);
	
	std::vector<KummerLineIsogenyGeneric*> psi_list;
	std::vector<KummerLineIsogenyGeneric*> phi_list;

	std::vector<std::array<fmpz_t*, 2>> factors = Helper::factoring(cofactor);

	for (size_t i = 0; i < factors.size(); ++i) {

		// Map P through chain lenght e of l-isogenies
		P = evaluate_factored_kummer_isogeny(psi_list, P);
		psi_list.clear();

		// compute point Q of order l^e
		fmpz_pow_fmpz(temp, *factors[i][0], *factors[i][1]);
		fmpz_divexact(cofactor, cofactor, temp);

		// use Q as kernel of degree l^e isogeny
		psi_list = sparse_isogeny_prime_power(P.multiply(cofactor), *factors[i][0], *factors[i][1], threshold);

		phi_list.insert(phi_list.begin(), psi_list.begin(), psi_list.end());
	}

	return phi_list;
}

KummerLineIsogenyComposite::KummerLineIsogenyComposite(KummerLine& domain, KummerPoint& kernel, fmpz_t& degree, int threshold) {
	// constructor for KummerLineIsogenyComposite class

	i_phis = EvaluateIsogenies::factored_kummer_isogeny(domain, kernel, degree, threshold);

	fmpz_set(i_degree, i_phis[0]->i_degree);
	for (size_t i = 1; i < i_phis.size(); ++i) {
		fmpz_mul(i_degree, i_degree, i_phis[i]->i_degree);
	}

	i_domain = i_phis[0]->i_domain;
	i_codomain = i_phis[i_phis.size() - 1]->i_codomain;
}

KummerPoint KummerLineIsogenyComposite::call(KummerPoint& P) { 
	// evaluates a KummerPoint through a chain of isogenies associated with the class

	return EvaluateIsogenies::evaluate_factored_kummer_isogeny(i_phis, P); 
}

KummerLineIsogeny::KummerLineIsogeny(KummerLine& domain, KummerPoint& kernel, fmpz_t degree) {
	// constructor for KummerLineIsogeny class

	fmpz_set(i_degree, degree);
	i_kernel = kernel;
	i_domain = domain;
	i_codomain = compute_codomain();
}

KummerPoint KummerLineIsogeny::call(KummerPoint& P) {
	// evaluates isogeny of a point P

	if (fmpz_cmp_ui(i_degree, 2) == 0) {
		return evaluate_isogeny_even(P);
	}
	return evaluate_isogeny(P);
}

std::vector<std::array<fq_poly_t*, 2>> KummerLineIsogeny::precompute_edwards_multiples(fmpz_t& d) {
	// compute [i]k for i in [1...d]

	fmpz_t i;
	std::vector<KummerPoint> K_muls = i_kernel.multiples();
	std::vector<std::array<fq_poly_t*, 2>> E_muls;
	for (fmpz_zero(i); fmpz_cmp(i, d) < 0; fmpz_add_ui(i, i, 1)) {
		fq_poly_t YE, ZE;
		fq_poly_sub(YE, K_muls[fmpz_get_ui(i)].p_x, K_muls[fmpz_get_ui(i)].p_z, *K_muls[fmpz_get_ui(i)].p_parent.field);
		fq_poly_add(ZE, K_muls[fmpz_get_ui(i)].p_x, K_muls[fmpz_get_ui(i)].p_z, *K_muls[fmpz_get_ui(i)].p_parent.field);
		E_muls.push_back({ &YE, &ZE });
	}

	// deallocate unused dynamic memory
	fmpz_clear(i);
	K_muls.clear();

	return E_muls;
}

std::array<fq_poly_t*, 2> KummerLineIsogeny::compute_codomain_constants() {
	// compute codomain constants

	// compute and store pairs of points for later evaluation
	fmpz_t d;
	fmpz_sub_ui(d, i_degree, 1);
	fmpz_divexact_ui(d, d, 2);
	i_edwards_multiples = precompute_edwards_multiples(d);

	// convert to twisted Edwards curve paramaters
	fq_poly_t Ded, Aed, A, C;

	fq_poly_set(A, i_domain.k_a, *i_domain.field);
	fq_poly_set(C, i_domain.k_c, *i_domain.field);
	fq_poly_add(Ded, C, C, *i_domain.field);
	fq_poly_add(Aed, A, Ded, *i_domain.field);
	fq_poly_sub(Ded, A, Ded, *i_domain.field);

	fq_poly_t prod_Y, prod_Z;
	fq_poly_one(prod_Y, *i_domain.field);
	fq_poly_one(prod_Z, *i_domain.field);
	for (size_t i = 0; i < i_edwards_multiples.size(); ++i) {
		fq_poly_mul(prod_Y, prod_Y, *i_edwards_multiples[i][0], *i_domain.field);
		fq_poly_mul(prod_Z, prod_Z, *i_edwards_multiples[i][1], *i_domain.field);
	}

	fq_poly_pow(prod_Y, prod_Y, 2, *i_domain.field);
	fq_poly_pow(prod_Y, prod_Y, 2, *i_domain.field);
	fq_poly_pow(prod_Y, prod_Y, 2, *i_domain.field);

	fq_poly_pow(prod_Z, prod_Z, 2, *i_domain.field);
	fq_poly_pow(prod_Z, prod_Z, 2, *i_domain.field);
	fq_poly_pow(prod_Z, prod_Z, 2, *i_domain.field);

	fq_poly_pow(Aed, Aed, fmpz_get_ui(i_degree), *i_domain.field);
	fq_poly_mul(Aed, Aed, prod_Z, *i_domain.field);

	fq_poly_pow(Ded, Ded, fmpz_get_ui(i_degree), *i_domain.field);
	fq_poly_mul(Ded, Ded, prod_Y, *i_domain.field);

	// change back to montgomery paramaters
	fq_poly_add(A, Aed, Ded, *i_domain.field);
	fq_poly_sub(C, Aed, Ded, *i_domain.field);
	fq_poly_add(A, A, A, *i_domain.field);

	// deallocate unused dynamic memory
	fmpz_clear(d);
	fq_poly_clear(Ded, *i_domain.field);
	fq_poly_clear(Aed, *i_domain.field);
	fq_poly_clear(A, *i_domain.field);
	fq_poly_clear(C, *i_domain.field);
	fq_poly_clear(prod_Y, *i_domain.field);
	fq_poly_clear(prod_Z, *i_domain.field);

	return { &A, &C };
}

std::array<fq_poly_t*, 2> KummerLineIsogeny::compute_codomain_constants_even() {
	// compute codomain constants of an even degree

	// C = ZK^2
	fq_poly_t A, C, XK, ZK;
	fq_poly_set(XK, i_kernel.p_x, *i_kernel.p_parent.field);
	fq_poly_set(ZK, i_kernel.p_z, *i_kernel.p_parent.field);
	fq_poly_mul(C, ZK, ZK, *i_kernel.p_parent.field);

	// A = 2*(ZK^2 - 2*XK^2)
	fq_poly_mul(A, XK, XK, *i_kernel.p_parent.field);
	fq_poly_add(A, A, A, *i_kernel.p_parent.field);
	fq_poly_sub(A, C, A, *i_kernel.p_parent.field);
	fq_poly_add(A, A, A, *i_kernel.p_parent.field);

	// deallocate unused dynamic memory
	fq_poly_clear(XK, *i_kernel.p_parent.field);
	fq_poly_clear(ZK, *i_kernel.p_parent.field);

	return { &A, &C };
}

KummerLine KummerLineIsogeny::compute_codomain() {
	// computes codomain constants

	std::array<fq_poly_t*, 2> codomains;
	if (fmpz_cmp_ui(i_degree, 2) == 0) {
		codomains = compute_codomain_constants_even();
	}
	else {
		codomains = compute_codomain_constants();
	}
	return KummerLine(*i_domain.field, codomains);
}

KummerPoint KummerLineIsogeny::evaluate_isogeny(KummerPoint& P) {
	// evaluates isogeny of point P

	fq_poly_t XP, ZP, Psum, Pdiff, X_new, Z_new;

	fq_poly_set(XP, P.p_x, *P.p_parent.field);
	fq_poly_set(ZP, P.p_z, *P.p_parent.field);
	fq_poly_add(Psum, XP, ZP, *P.p_parent.field);
	fq_poly_sub(Pdiff, XP, ZP, *P.p_parent.field);

	// loop through the d-multiples
	fq_poly_one(X_new, *P.p_parent.field);
	fq_poly_one(Z_new, *P.p_parent.field);
	
	fq_poly_t diff_EZ, sum_EY, temp1, temp2;
	for (size_t i = 0; i < i_edwards_multiples.size(); ++i) {
		fq_poly_mul(diff_EZ, Pdiff, *i_edwards_multiples[i][1], *P.p_parent.field);
		fq_poly_mul(sum_EY, *i_edwards_multiples[i][0], Psum, *P.p_parent.field);
		fq_poly_add(temp1, diff_EZ, sum_EY, *P.p_parent.field);
		fq_poly_sub(temp2, diff_EZ, sum_EY, *P.p_parent.field);
		fq_poly_mul(X_new, X_new, temp1, *P.p_parent.field);
		fq_poly_mul(Z_new, Z_new, temp2, *P.p_parent.field);
	}

	// square and multiply with original
	fq_poly_pow(X_new, X_new, 2, *P.p_parent.field);
	fq_poly_mul(X_new, X_new, XP, *P.p_parent.field);
	fq_poly_pow(Z_new, Z_new, 2, *P.p_parent.field);
	fq_poly_mul(Z_new, Z_new, ZP, *P.p_parent.field);

	// deallocate unused dynamic memory
	fq_poly_clear(XP, *P.p_parent.field);
	fq_poly_clear(ZP, *P.p_parent.field);
	fq_poly_clear(Psum, *P.p_parent.field);
	fq_poly_clear(Pdiff, *P.p_parent.field);
	fq_poly_clear(diff_EZ, *P.p_parent.field);
	fq_poly_clear(sum_EY, *P.p_parent.field);
	fq_poly_clear(temp1, *P.p_parent.field);
	fq_poly_clear(temp2, *P.p_parent.field);

	return KummerPoint(i_codomain, { &X_new, &Z_new });
}

KummerPoint KummerLineIsogeny::evaluate_isogeny_even(KummerPoint& P) {
	// evaluates isogeny of point P of even degree

	fq_poly_t XK, ZK, XP, ZP, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9;

	fq_poly_set(XK, i_kernel.p_x, *i_kernel.p_parent.field);
	fq_poly_set(ZK, i_kernel.p_z, *i_kernel.p_parent.field);
	fq_poly_set(XP, P.p_x, *P.p_parent.field);
	fq_poly_set(ZP, P.p_x, *P.p_parent.field);

	fq_poly_add(T0, XK, ZK, *P.p_parent.field);
	fq_poly_sub(T1, XK, ZK, *P.p_parent.field);
	fq_poly_add(T2, XK, ZP, *P.p_parent.field);
	fq_poly_sub(T3, ZP, XP, *P.p_parent.field);
	fq_poly_mul(T4, T3, T0, *P.p_parent.field);
	fq_poly_mul(T5, T2, T1, *P.p_parent.field);
	fq_poly_sub(T6, T4, T5, *P.p_parent.field);
	fq_poly_add(T7, T4, T5, *P.p_parent.field);
	fq_poly_mul(T8, XP, T6, *P.p_parent.field); // XP * ((ZP - XP)(XK + ZK) - (XP + ZP)(XK - ZK))
	fq_poly_mul(T9, ZP, T7, *P.p_parent.field); // ZP * ((ZP - XP)(XK + ZK) + (XP + ZP)(XK - ZK))

	// deallocate unused dynamic memory
	fq_poly_clear(XK, *P.p_parent.field);
	fq_poly_clear(ZK, *P.p_parent.field);
	fq_poly_clear(XP, *P.p_parent.field);
	fq_poly_clear(ZP, *P.p_parent.field);
	fq_poly_clear(T0, *P.p_parent.field);
	fq_poly_clear(T1, *P.p_parent.field);
	fq_poly_clear(T2, *P.p_parent.field);
	fq_poly_clear(T3, *P.p_parent.field);
	fq_poly_clear(T4, *P.p_parent.field);
	fq_poly_clear(T5, *P.p_parent.field);
	fq_poly_clear(T6, *P.p_parent.field);
	fq_poly_clear(T7, *P.p_parent.field);

	return KummerPoint(i_codomain, { &T8, &T9 });
}

KummerLineIsogeny_VeluSqrt::KummerLineIsogeny_VeluSqrt(KummerLine& domain, KummerPoint& kernel, fmpz_t degree) {
	// constructor for KummerLineIsogeny_VeluSqrt

	fmpz_set(i_degree, degree);
	i_kernel = kernel;
	i_domain = domain;
	fq_poly_set(i_a, *domain.a(), *domain.field);
	i_ring = domain.field;
	fq_poly_gen(i_ring_generator, *i_ring);
	fmpz_t b, c, temp;
	fmpz_sub_ui(temp, degree, 1);
	fmpz_sqrt(b, temp);
	fmpz_divexact_ui(b, b, 2);
	fmpz_mul_ui(c, b, 4);
	fmpz_divexact(c, temp, c);
	fmpz_mul_ui(b, b, 4);
	fmpz_mul(b, b, c);
	fmpz_sub(stop, degree, b);
	
	hI_tree = hI_precomputation(b, c);
	EJ_parts = EJ_precomputation(b);
	hK = hK_precomputation(b, c);

	i_codomain = compute_codomain();

	// deallocate unused dynamic memory
	fmpz_clear(degree);
	fmpz_clear(b);
	fmpz_clear(c);
	fmpz_clear(temp);
}

fmpz_t* KummerLineIsogeny_VeluSqrt::hI_resultant(fq_poly_t& poly) {
	// calculates the polynomial resultant from a polynomial product tree and a polynomial

	return ProductTree::product_tree_resultant(hI_tree, poly);
}

ProductTree KummerLineIsogeny_VeluSqrt::hI_precomputation(fmpz_t& b, fmpz_t& c) {
	// precompute h_I using product tree

	fmpz_t temp;
	fmpz_mul_ui(temp, b, 2);
	KummerPoint Q = i_kernel.multiply(temp);
	KummerPoint step = Q.twice();
	KummerPoint diff = Q;
	KummerPoint temp1;
	KummerPoint temp2;
	std::vector<fq_poly_t*> leaves;
	for (fmpz_zero(temp); fmpz_cmp(temp, c) < 0; fmpz_add_ui(temp, temp, 1)) {
		fq_poly_t next;
		fq_poly_sub(next, i_ring_generator, *Q.x(), *i_ring);
		leaves.push_back(&next);
		if (fmpz_cmp(temp, c-1) < 0) {
			temp1 = Q.add(step, diff);
			temp2 = Q;
			Q = temp1;
			diff = temp2;
		}
	}

	// deallocate dynamic memory
	fmpz_clear(temp);
	delete& Q;
	delete& step;
	delete& diff;
	delete& temp1;
	delete& temp2;


	return ProductTree(leaves, i_ring);
}

std::array<fq_poly_t*, 3> KummerLineIsogeny_VeluSqrt::Fs(fq_poly_t& X1, fq_poly_t& X2) {
	// calculates the biquadratic polynomials associated with X1 and X2

	fq_poly_t X1X2, param1, param2, param3, temp;
	fq_t multiplier;
	fq_poly_mul(X1X2, X1, X2, *i_ring);

	// param1 computation
	fq_poly_sub(param1, X1, X2, *i_ring);
	fq_poly_pow(param1, param1, 2, *i_ring);

	// param2 computation
	fq_poly_add_si(param2, X1X2, 1, *i_ring);
	fq_poly_add(temp, X1, X2, *i_ring);
	fq_poly_mul(param2, param2, temp, *i_ring);
	fq_poly_mul(temp, i_a, X1X2, *i_ring);
	fq_set_si(multiplier, 2, *i_ring);
	fq_poly_scalar_mul_fq(temp, temp, multiplier, *i_ring);
	fq_poly_add(param2, param2, temp, *i_ring);
	fq_set_si(multiplier, -2, *i_ring);
	fq_poly_scalar_mul_fq(param2, param2, multiplier, *i_ring);

	// param3 computation
	fq_poly_add_si(param3, X1X2, -1, *i_ring);
	fq_poly_pow(param3, param3, 2, *i_ring);

	std::array<fq_poly_t*, 3> polys = { &param1, &param2, &param3 };
	 
	// deallocate unused dynamic memory
	fq_clear(multiplier, *i_ring);
	fq_poly_clear(X1X2, *i_ring);
	fq_poly_clear(temp, *i_ring);

	return polys;
}

std::vector<std::array<fq_poly_t*, 3>> KummerLineIsogeny_VeluSqrt::EJ_precomputation(fmpz_t& b) {
	// precompute EJ parts

	KummerPoint Q = KummerPoint(i_kernel.p_parent, { &i_kernel.p_x, &i_kernel.p_z });
	KummerPoint step = Q.twice();
	KummerPoint diff = Q;
	KummerPoint temp1;
	KummerPoint temp2;
	fmpz_t i;
	std::vector<std::array<fq_poly_t*, 3>> EJ_parts;
	for (fmpz_zero(i); fmpz_cmp(i, b) < 0; fmpz_add_ui(i, i, 1)) {
		EJ_parts.push_back(Fs(i_ring_generator, *Q.x()));
		if (i < b - 1){
			temp1 = Q.add(step, diff);
			temp2 = Q;
			Q = temp1;
			diff = temp2;
		}
	}

	// deallocate unused dynamic memory
	fmpz_clear(i);
	delete& Q;
	delete& step;
	delete& diff;
	delete& temp1;
	delete& temp2;
	return EJ_parts;
}

fq_poly_t* KummerLineIsogeny_VeluSqrt::hK_precomputation(fmpz_t& b, fmpz_t& c) {
	// precompute hk

	std::vector<fq_poly_t*> hK;
	KummerPoint Q = i_kernel.twice();
	KummerPoint step = Q;
	KummerPoint next_point = Q.twice();
	KummerPoint temp1;
	KummerPoint temp2;
	fmpz_t stop, i;
	fmpz_mul(stop, b, c);
	fmpz_mul_ui(stop, stop, 4);
	fmpz_sub(stop, i_degree, stop);
	for (fmpz_set_ui(i, 2); fmpz_cmp(i, stop) < 0; fmpz_add_ui(i, i, 2)) {
		fq_poly_t each;
		fq_poly_mul(each, Q.p_z, i_ring_generator, *i_ring);
		fq_poly_sub(each, each, Q.p_x, *i_ring);
		hK.push_back(&each);
		if (i < stop - 1) {
			temp1 = next_point;
			temp2 = next_point.add(step, Q);
			Q = temp1;
			next_point = temp2;
		}
	}

	fq_poly_t* result = ProductTree::prod(hK, i_ring);

	// deallocate unused dynamic memory
	fmpz_clear(stop);
	fmpz_clear(i);
	hK.clear();
	delete& Q;
	delete& step;
	delete& next_point;
	delete& temp1;
	delete& temp2;

	return result;
}

std::array<fq_poly_t*, 2> KummerLineIsogeny_VeluSqrt::compute_codomain_constants() {
	// computes codomain constants

	// polynomials for alpha = 1 and alpha = -1
	std::vector<fq_poly_t*> values1;
	std::vector<fq_poly_t*> values2;
	for (size_t i = 0; i < EJ_parts.size(); ++i) {
		fq_poly_t E0, E1;
		fq_poly_add(E0, *EJ_parts[i][0], *EJ_parts[i][1], *i_ring);
		fq_poly_add(E0, E0, *EJ_parts[i][2], *i_ring);
		fq_poly_sub(E1, *EJ_parts[i][0], *EJ_parts[i][1], *i_ring);
		fq_poly_add(E1, E1, *EJ_parts[i][2], *i_ring);
		values1.push_back(&E0);
		values2.push_back(&E1);
	}
	fq_poly_t* E0J = ProductTree::prod(values1, i_ring);
	fq_poly_t* E1J = ProductTree::prod(values2, i_ring);

	// compute resultants and evaluate hK at 1 and -1
	fmpz_t* R0 = hI_resultant(*E0J);
	fmpz_t* R1 = hI_resultant(*E1J);

	// compute num and den
	fq_t M0, M1, num, den;
	fq_one(M0, *i_ring);
	fq_set_si(M1, -1, *i_ring);
	fq_poly_evaluate_fq(M0, *hK, M0, *i_ring);
	fq_poly_evaluate_fq(M1, *hK, M1, *i_ring);

	fq_mul_fmpz(num, M0, *R0, *i_ring);
	fq_mul_fmpz(den, M1, *R1, *i_ring);

	fmpz_t power;
	fmpz_set_ui(power, 2);
	for (int i = 0; i < 3; ++i) {
		fq_pow(num, num, power, *i_ring);
		fq_pow(den, den, power, *i_ring);
	}

	fq_poly_t numerator, denominator;

	// calculate numerator
	fq_poly_add_si(numerator, i_a, -2, *i_ring);
	fq_poly_pow(numerator, numerator, *i_degree, *i_ring);
	fq_poly_scalar_mul_fq(numerator, numerator, num, *i_ring);

	//calculate denominator
	fq_poly_add_si(denominator, i_a, 2, *i_ring);
	fq_poly_pow(denominator, denominator, *i_degree, *i_ring);
	fq_poly_scalar_mul_fq(denominator, denominator, den, *i_ring);

	// calculate the new curve constants
	fq_poly_t A_new, C_new;

	fq_poly_add(A_new, numerator, denominator, *i_ring);
	fq_poly_sub(C_new, denominator, numerator, *i_ring);
	fq_poly_add(A_new, A_new, A_new, *i_ring);

	// deallocate unused dynamic memory
	fq_poly_clear(*E0J, *i_ring);
	fq_poly_clear(*E1J, *i_ring);
	values1.clear();
	values2.clear();
	fmpz_clear(*R0);
	fmpz_clear(*R1);
	fq_clear(M0, *i_ring);
	fq_clear(M1, *i_ring);
	fq_clear(num, *i_ring);
	fq_clear(den, *i_ring);
	fmpz_clear(power);
	fq_poly_clear(numerator, *i_ring);
	fq_poly_clear(denominator, *i_ring);

	return { &A_new, &C_new };
}

KummerLine KummerLineIsogeny_VeluSqrt::compute_codomain() {
	// computes codomain constants

	std::array<fq_poly_t*, 2> AC = compute_codomain_constants();
	return KummerLine(*i_ring, AC);
}

KummerPoint KummerLineIsogeny_VeluSqrt::evaluate_isogeny(KummerPoint& P) {
	// evaluates isogeny of a Point P

	// checks for point at infinity
	if (fq_poly_is_zero(P.p_z, *P.p_parent.field)) {
		fq_poly_t one, zero;
		fq_poly_one(one, *P.p_parent.field);
		fq_poly_one(zero, *P.p_parent.field);
		return KummerPoint(i_codomain, { &one, &zero });
	}

	// x coordinate of point to evaluate
	fq_poly_t* alpha = P.x();
	fq_poly_t alpha_inv, r, one;
	fq_poly_one(one, *i_ring);
	fq_poly_divrem(alpha_inv, r, one, *alpha, *i_ring);

	std::vector<fq_poly_t*> values1;
	std::vector<fq_poly_t*> values2;

	for (size_t i = 0; i < EJ_parts.size(); ++i) {
		fq_poly_t E0, E1;
		fq_poly_mul(E0, *EJ_parts[i][0], alpha_inv, *i_ring);
		fq_poly_add(E0, E0, *EJ_parts[i][1], *i_ring);
		fq_poly_mul(E0, E0, alpha_inv, *i_ring);
		fq_poly_add(E0, E0, *EJ_parts[i][2], *i_ring);

		fq_poly_mul(E1, *EJ_parts[i][0], *alpha, *i_ring);
		fq_poly_add(E1, E1, *EJ_parts[i][1], *i_ring);
		fq_poly_mul(E1, E1, *alpha, *i_ring);
		fq_poly_add(E1, E1, *EJ_parts[i][2], *i_ring);
		values1.push_back(&E0);
		values2.push_back(&E1);
	}

	// calculate resultants
	fq_poly_t* EJ0 = ProductTree::prod(values1, i_ring);
	fq_poly_t* EJ1 = ProductTree::prod(values2, i_ring);

	fmpz_t* R0 = hI_resultant(*EJ0);
	fmpz_t* R1 = hI_resultant(*EJ1);


	// calculate evalutations
	fq_t coeff1, coeff0;
	fq_poly_t M0, M1;

	fq_poly_get_coeff(coeff1, *hK, 1, *i_ring);
	fq_poly_get_coeff(coeff0, *hK, 0, *i_ring);
	fq_poly_scalar_mul_fq(M0, alpha_inv, coeff1, *i_ring);
	fq_poly_scalar_addmul_fq(M0, one, coeff0, *i_ring);

	fq_poly_get_coeff(coeff1, *hK, 1, *i_ring);
	fq_poly_get_coeff(coeff0, *hK, 0, *i_ring);
	fq_poly_scalar_mul_fq(M1, *alpha, coeff1, *i_ring);
	fq_poly_scalar_addmul_fq(M1, one, coeff0, *i_ring);

	// make new point
	fq_poly_t X_new, Z_new;

	fq_set_fmpz(coeff0, *R0, *i_ring);
	fq_set_fmpz(coeff1, *R1, *i_ring);
	fq_poly_pow(Z_new, *alpha, *i_degree, *i_ring);
	fq_poly_scalar_mul_fq(X_new, M0, coeff0, *i_ring);
	fq_poly_pow(X_new, X_new, 2, *i_ring);
	fq_poly_mul(X_new, X_new, Z_new, *i_ring);
	fq_poly_scalar_mul_fq(Z_new, M1, coeff1, *i_ring);
	fq_poly_pow(Z_new, Z_new, 2, *i_ring);

	// deallocate unused dynamic memory
	fq_poly_clear(*alpha, *i_ring);
	fq_poly_clear(alpha_inv, *i_ring);
	fq_poly_clear(r, *i_ring);
	fq_poly_clear(one, *i_ring);
	values1.clear();
	values2.clear();
	fq_poly_clear(*EJ0, *i_ring);
	fq_poly_clear(*EJ1, *i_ring);
	fmpz_clear(*R0);
	fmpz_clear(*R1);
	fq_clear(coeff1, *i_ring);
	fq_clear(coeff0, *i_ring);
	fq_poly_clear(M0, *i_ring);
	fq_poly_clear(M1, *i_ring);

	return KummerPoint(i_codomain, { &X_new, &Z_new });
}