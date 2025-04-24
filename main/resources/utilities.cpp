#include "../headers/utilities.h"

void TorsionData::batch_cofactor_mul_generic(std::vector<Point>& G_list, fq_ctx_t& field, std::vector<slong>& pis, fmpz_t lower, fmpz_t upper) {
	// applies the scalar multiplication to elements in G_list from prime factors of group order

	// last recursion step does not need further splitting
	fmpz_t test;
	fmpz_sub(test, upper, lower);
	if (fmpz_is_one(test)) {
		fmpz_clear(lower);
		fmpz_clear(upper);
		fmpz_clear(test);
		return;
	}

	// split list into two parts
	fmpz_t mid, cl, cu, i;
	fmpz_sub(mid, upper, lower);
	fmpz_add_ui(mid, mid, 1);
	fmpz_divexact_ui(mid, mid, 2);
	fmpz_add(mid, mid, lower);
	fmpz_one(cl);
	fmpz_one(cu);
	for (fmpz_set(i, lower); fmpz_cmp(i, mid) < 0; fmpz_add_ui(i, i, 1)) {
		fmpz_mul_ui(cu, cu, pis[fmpz_get_ui(i)]);
	}
	for (fmpz_set(i, mid); fmpz_cmp(i, upper) < 0; fmpz_add_ui(i, i, 1)) {
		fmpz_mul_ui(cl, cl, pis[fmpz_get_ui(i)]);
	}

	// multiply to get new start points for two sublists
	G_list[fmpz_get_ui(mid)] = Point::scalar_multiplication(G_list[fmpz_get_ui(lower)], cu);
	G_list[fmpz_get_ui(lower)] = Point::scalar_multiplication(G_list[fmpz_get_ui(lower)], cl);

	// call function recursively for each sublist
	batch_cofactor_mul_generic(G_list, field, pis, lower, mid);
	batch_cofactor_mul_generic(G_list, field, pis, mid, upper);

	// deallocate dynamic memory
	fmpz_clear(mid);
	fmpz_clear(cl);
	fmpz_clear(cu);
	fmpz_clear(i);
	fmpz_clear(lower);
	fmpz_clear(upper);
	fmpz_clear(test);
}

void TorsionData::batch_cofactor_mul_generic(std::vector<fq_poly_t*>& G_list, fq_ctx_t& field, std::vector<slong>& pis, fmpz_t lower, fmpz_t upper) {
	// performs polynomials to the power of prime factors of group order

	// last recursion step does not require further splitting
	fmpz_t test;
	fmpz_sub(test, upper, lower);
	if (fmpz_is_one(test)) {
		fmpz_clear(lower);
		fmpz_clear(upper);
		fmpz_clear(test);
		return;
	}

	// split list into two parts
	fmpz_t mid, cl, cu, i;
	fmpz_sub(mid, upper, lower);
	fmpz_add_ui(mid, mid, 1);
	fmpz_divexact_ui(mid, mid, 2);
	fmpz_add(mid, mid, lower);
	fmpz_one(cl);
	fmpz_one(cu);
	for (fmpz_set(i, lower); fmpz_cmp(i, mid) < 0; fmpz_add_ui(i, i, 1)) {
		fmpz_mul_ui(cu, cu, pis[fmpz_get_ui(i)]);
	}
	for (fmpz_set(i, mid); fmpz_cmp(i, upper) < 0; fmpz_add_ui(i, i, 1)) {
		fmpz_mul_ui(cl, cl, pis[fmpz_get_ui(i)]);
	}

	// multiply to get new start points for two sublists
	fq_poly_pow(*G_list[fmpz_get_ui(mid)], *G_list[fmpz_get_ui(mid)], fmpz_get_ui(cu), field);
	fq_poly_pow(*G_list[fmpz_get_ui(lower)], *G_list[fmpz_get_ui(lower)], fmpz_get_ui(cl), field);

	// call function recursively for each sublist
	batch_cofactor_mul_generic(G_list, field, pis, lower, mid);
	batch_cofactor_mul_generic(G_list, field, pis, mid, upper);

	// deallocate dynamic memory
	fmpz_clear(mid);
	fmpz_clear(cl);
	fmpz_clear(cu);
	fmpz_clear(i);
	fmpz_clear(lower);
	fmpz_clear(upper);
	fmpz_clear(test);
}

std::vector<slong> TorsionData::has_order_constants(fmpz_t& D) {
	// helper function, finds constants to help with has_order_D

	std::vector<std::array<fmpz_t*, 2>> factors = Helper::factoring(D);
	std::vector<slong> pis;
	for (size_t i = 0; i < factors.size(); ++i) {
		for (slong x = 0; x < fmpz_get_si(*factors[i][1]); ++x) {
			pis.push_back(fmpz_get_si(*factors[i][0]));
		}
	}
	fmpz_t* D_radical = ProductTree::prod(pis);
	fmpz_divexact(D, D, *D_radical);
	factors.clear();
	return pis;
}

bool TorsionData::has_order_D(Point G, fmpz_t& D) {
	// checks a point G has group order D

	// checks G isn't additive identity 
	std::cout << "has order D";
	if (fq_poly_is_zero(G.p_z, *G.curve.field)) {
		return false;
	}

	fmpz_t D_top;
	fmpz_set(D_top, D);

	std::vector<slong> pis = has_order_constants(D_top);

	// abort early if Gtop reveals identity
	Point Gtop = Point::scalar_multiplication(G, D_top);
	if (fq_poly_is_zero(Gtop.p_z, *G.curve.field)) {
		return false;
	}

	// create G_list, where the first element is Gtop and the rest are at points x=0, y=0
	std::vector<Point> G_list;
	G_list.push_back(Gtop);

	fq_poly_t zero;
	fq_poly_t one;

	fq_poly_zero(zero, *G.curve.field);
	fq_poly_one(one, *G.curve.field);

	for (size_t i = 0; i < pis.size(); ++i) {
		G_list.push_back(Point(zero, zero, one, G.curve));
	}

	// checks identity can't be revealed by removing factors of the order of the group

	fmpz_t start;
	fmpz_t end;
	fmpz_zero(start);
	fmpz_set_ui(end, pis.size());
	if (pis.size() > 1) {
		batch_cofactor_mul_generic(G_list, *G.curve.field, pis, start, end);
		for (size_t i = 0; i < G_list.size(); ++i) {
			if (fq_poly_is_zero(G_list[i].p_z, *G.curve.field)) {
				fmpz_clear(D_top);
				pis.clear();
				delete &Gtop;
				G_list.clear();
				fq_poly_clear(one, *G.curve.field);
				fq_poly_clear(zero, *G.curve.field);
				fmpz_clear(start);
				fmpz_clear(end);
				return false;
			}
		}
	}

	// memory cleanup
	fmpz_clear(D_top);
	pis.clear();
	delete& Gtop;
	G_list.clear();
	fq_poly_clear(one, *G.curve.field);
	fq_poly_clear(zero, *G.curve.field);
	fmpz_clear(start);
	fmpz_clear(end);
	return true;
}

bool TorsionData::has_order_D(fq_poly_t& G, fq_ctx_t& field, fmpz_t& D) {
	// checks a group element G over field has order D

	// checks G isn't multiplicative identity element
	if (fq_poly_is_one(G, field)) {
		return false;
	}

	fmpz_t D_top;
	fmpz_set(D_top, D);

	std::vector<slong> pis = has_order_constants(D_top);

	// abort early if Gtop reveals identity
	fq_poly_t Gtop;
	fq_poly_pow(Gtop, G, fmpz_get_si(D_top), field);
	if (fq_poly_is_one(G, field)) {
		fmpz_clear(D_top);
		fq_poly_clear(Gtop, field);
		return false;
	}

	// create G_list where Gtop is G and the rest are polynomial 1
	std::vector<fq_poly_t*> G_list;
	G_list.push_back(&Gtop);

	for (size_t i = 0; i < pis.size() - 1; ++i) {
		fq_poly_t next;
		fq_poly_one(next, field);
		G_list.push_back(&next);
	}

	// check identity can't be revealed from removing the factors of the order of the group
	fmpz_t start;
	fmpz_t end;
	fmpz_zero(start);
	fmpz_set_ui(end, pis.size());

	if (pis.size() > 1) {
		batch_cofactor_mul_generic(G_list, field, pis, start, end);
		for (size_t i = 0; i < G_list.size(); ++i) {
			if (fq_poly_is_one(*G_list[i], field)) {
				fmpz_clear(D_top);
				fq_poly_clear(Gtop, field);
				pis.clear();
				G_list.clear();
				fmpz_clear(start);
				fmpz_clear(end);
				return false;
			}
		}
	}
	// cleanup memory
	fmpz_clear(D_top);
	fq_poly_clear(Gtop, field);
	pis.clear();
	G_list.clear();
	fmpz_clear(start);
	fmpz_clear(end);
	return true;
}

Point TorsionData::generate_random_point(EllipticCurve& E, flint_rand_t& prng){
	// generates random point on an elliptic curve

	std::cout << "generating point";
	Point P = Point::generate_random_element(E, prng);
	std::cout << "generated point";
	fq_poly_t neg_y;
	fq_poly_init(neg_y, *E.field);
	slong degree1 = fq_poly_degree(P.p_y, *E.field);
	slong degree2 = fq_poly_degree(neg_y, *E.field);

	// finds the negative version of point P
	fq_poly_neg(neg_y, P.p_y, *P.curve.field);

	// compares and returns the point with the smallest y coordinate
	if (degree1 > degree2){
		fq_poly_set(P.p_y, neg_y, *E.field);
		fq_poly_clear(neg_y, *E.field);
	}
	else if (degree1 < degree2) {
		fq_poly_clear(neg_y, *E.field);
	}
	else{
		fq_t coeff1, coeff2;
		fmpz_t cmp1, cmp2;
		// compare y coefficients from highest degree to lowest and returns Point with the smallest y
		for (slong i = degree1; i >= 0; --i) {

			fq_poly_get_coeff(coeff1, P.p_y, i, *E.field);
			fq_poly_get_coeff(coeff2, neg_y, i, *E.field);

			fq_get_fmpz(cmp1, coeff1, *E.field);
			fq_get_fmpz(cmp2, coeff1, *E.field);

			int cmp = fmpz_cmp(cmp1, cmp2);
			if (cmp < 0) {
				fq_clear(coeff1, *E.field);
				fq_clear(coeff2, *E.field);
				fq_poly_clear(neg_y, *E.field);
				fmpz_clear(cmp1);
				fmpz_clear(cmp2);
				return P;
			}
			else if (cmp > 0) {
				fq_poly_set(P.p_y, neg_y, *E.field);
				fq_clear(coeff1, *E.field);
				fq_clear(coeff2, *E.field);
				fq_poly_clear(neg_y, *E.field);
				fmpz_clear(cmp1);
				fmpz_clear(cmp2);
				return P;
			}
		}
		fq_clear(coeff1, *E.field);
		fq_clear(coeff2, *E.field);
		fq_poly_clear(neg_y, *E.field);
		fmpz_clear(cmp1);
		fmpz_clear(cmp2);
	}
	return P;
}

Point TorsionData::generate_point_order_D(EllipticCurve& E, fmpz_t& D, fmpz_t& p, flint_rand_t& prng){
	// generates a point of order D

	fmpz_t n;
	fmpz_add_ui(n, p, 1);
	fmpz_divexact(n, n, D);

	for (int i = 0; i < 1000; ++i) {
		Point P = generate_random_point(E, prng);
		P = Point::scalar_multiplication(P, n);

		// in case we randomly pick a point on the n-torsion
		if (fq_poly_is_zero(P.p_x, *E.field)) {
			continue;
		}

		// check that P has exactly order D
		if (has_order_D(P, D)) {
			fmpz_set(P.order, D);
			return P;
		}
	}
	throw std::runtime_error("Never found a point of order D");
}

Point TorsionData::generate_linearly_independant(EllipticCurve& E, Point& P, fmpz_t& D, fmpz_t& p, flint_rand_t& prng){
	// generates a point that is linearly independant from P

	for (int i = 0; i < 2000; ++i) {
		// generates a random point of order D
		Point Q = generate_point_order_D(E, D, p, prng);
		
		// checks to make sure the point is linearly independant
		fq_poly_t* pair = Point::weil_pairing(E, P, Q, D, prng);
		if (has_order_D(*pair, *P.curve.field, D)) {
			fmpz_set(Q.order, D);
			return Q;
		}
	}
	throw std::runtime_error("Never found a linearly independent point...");
}

std::array<Point, 2> TorsionData::torsion_basis(EllipticCurve& E, fmpz_t& D, fmpz_t& characteristic, flint_rand_t& prng) {
	// calculates the torsion points for the protocol

	// ensures D divides the curves order
	fmpz_t temp;
	fmpz_add_ui(temp, characteristic, 1);
	if (!fmpz_divisible(temp, D)){
		std::cout << "fail";
		throw std::runtime_error("D must divide the point's order");
	}

	// generate torsion data
	Point P = generate_point_order_D(E, D, characteristic, prng);
	Point Q = generate_linearly_independant(E, P, D, characteristic, prng);
	fmpz_clear(temp);
	return { P, Q };
}