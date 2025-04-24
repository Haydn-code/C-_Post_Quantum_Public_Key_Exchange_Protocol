#ifndef BIN_TERSIDH__HYBRID
#define BIN_TERSIDH__HYBRID

#include "bin-terSIDH.h"

class bin_terSIDH__hybrid {
	// class to run the hybrid versions of the bin and ter SIDH protocols
public:
	static std::array<fmpz_t, 4> compute_kernel_scalars_hybrid(std::vector<int> s);
	static std::tuple<fmpz_t*, KummerLine, KummerPoint, KummerPoint> keygenA(flint_rand_t& prng, fmpz_t& A, fmpz_t& B, KummerPoint& xQA, KummerPoint& xPA, KummerPoint& xPQA, KummerPoint& xPB, KummerPoint& xQB, KummerLine& L0);
	static std::tuple<std::array<fmpz_t, 4>, KummerLine, KummerPoint, KummerPoint, KummerPoint> keygenB(std::vector<int> sk_choices, int t, KummerPoint& xPB, KummerPoint& xQB, KummerLine& L0, KummerPoint& xPA, KummerPoint& xQA, KummerPoint& xPQA);
	static fq_poly_t* sharedA(fmpz_t& skA, KummerLine& EB, KummerPoint& RA, KummerPoint& SA, KummerPoint& RSA, fmpz_t& A);
	static fq_poly_t* sharedB(std::array<fmpz_t, 4> skB, KummerLine& EA, KummerPoint& RB, KummerPoint& SB);
	static void bin_terSIDH__hybrid_core(int t, std::vector<int> sk_choices, int l);
};

#endif