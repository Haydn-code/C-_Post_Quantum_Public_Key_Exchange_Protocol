#include "../headers/product_tree.h"

ProductTree::ProductTree(std::vector<fq_poly_t*> polys, fq_ctx_t* ring) {
	// constructor for product tree

	field = ring;

	layers.push_back(polys);

	// compute product tree layers
	while (layers.back().size() > 1) {
		std::vector<fq_poly_t*> next_layer;
		for (size_t i = 0; i < layers.back().size(); i += 2) {
			fq_poly_t temp;
			fq_poly_mul(temp, *layers.back()[i], *layers.back()[i], *field);
			next_layer.push_back(&temp);
		}
		layers.push_back(next_layer);
		next_layer.clear();
	}
}

std::vector<fq_poly_t*> ProductTree::remainders(fq_poly_t& x) {
	// function to compute remainders from a product tree

	std::vector<fq_poly_t*> X = { &x };

	fq_poly_t q;

	// loop through to find remainders
	for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
		std::vector<fq_poly_t*>& V = *it;

		std::vector<fq_poly_t*> nextX;
		for (size_t i = 0; i < V.size(); ++i) {
			fq_poly_t r;
			fq_poly_divrem(q, r, *X[i / 2], *V[i], *field);
			nextX.push_back(&r);
		}

		X = nextX;
		nextX.clear();
	}
	
	// deallocate dynamic memory
	fq_poly_clear(q, *field);

	return X;
}

fq_poly_t* ProductTree::prod(std::vector<fq_poly_t*>& polys, fq_ctx_t* field) {
	// function to calculate the product of a vector of polynomials over a finite field

	// set result to first item in polys
	fq_poly_t result;
	fq_poly_set(result, *polys[0], *field);

	// loop thorugh and multiply result by every poly in polys
	for (size_t i = 1; i < polys.size(); ++i) {
		fq_poly_mul(result, result, *polys[i], *field);
	}
	return &result;
}

fmpz_t* ProductTree::prod(std::vector<slong>& integers) {
	// function to calculate the product of a vector of integers

	// set result to the first item in integers
	fmpz_t result;
	fmpz_set_si(result, integers[0]);

	// loop through and multiply result by every integer in integers
	for (size_t i = 1; i < integers.size(); ++i) {
		fmpz_mul_si(result, result, integers[i]);
	}
	return &result;
}

fmpz_t* ProductTree::product_tree_resultant(ProductTree& hI_tree, fq_poly_t& poly) {
	// function to evaluate a resultant with `h_I` quickly using a product tree

	std::vector<fq_poly_t*> rems = hI_tree.remainders(poly);
	fq_poly_t* r = prod(rems, hI_tree.field);
	int s = 1;
	if (hI_tree.layers[0].size() % 2 == 1 and fq_poly_degree(*r, *hI_tree.field) == 1) {
		int s = -1;
	}
	fq_t x;
	fmpz_t result;
	fq_poly_get_coeff(x, *r, 0, *hI_tree.field);
	fq_get_fmpz(result, x, *hI_tree.field);
	fmpz_mul_ui(result, result, s);

	// deallocate dynamic memory
	fq_poly_clear(*r, *hI_tree.field);
	fq_clear(x, *hI_tree.field);

	return &result;
}