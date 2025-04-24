#ifndef PRODUCT_TREE_H
#define PRODUCT_TREE_H

#include <vector>
#include <tuple>
#include <flint/fmpz.h>
#include <flint/fq_poly.h>

class ProductTree {
	// Class for a binary product tree
public:
	std::vector<std::vector<fq_poly_t*>> layers;
	fq_ctx_t* field;
	ProductTree() { fmpz_t degree; fmpz_set_ui(degree, 2); fq_ctx_init(*field, degree, 2, "x"); fmpz_clear(degree); };
	ProductTree(std::vector<fq_poly_t*> polys, fq_ctx_t* ring);
	std::vector<fq_poly_t*> remainders(fq_poly_t& x);
	static fq_poly_t* prod(std::vector<fq_poly_t*>& polys, fq_ctx_t* field);
	static fmpz_t* prod(std::vector<slong>& integers);
	static fmpz_t* product_tree_resultant(ProductTree& hI_tree, fq_poly_t& poly);
};
#endif