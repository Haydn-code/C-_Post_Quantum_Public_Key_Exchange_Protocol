#include "headers/bin-terSIDH--hybrid.h"

std::array<fmpz_t, 4> bin_terSIDH__hybrid::compute_kernel_scalars_hybrid(std::vector<int> s) {
    /*
    given a ternary secret `s`, returns scalars `B0` and `B1`
    such that the isogeny associated with `s` and orientation (P, Q)
    has kernel <[B0]*P + [B1]*Q>.
    The function also returns `order0` and `order1`, the orders
    of points [B0]*P and [B1]*Q, which is used in the isogeny computations.
    */

    fmpz_t B0, B1, order0, order1;
    fmpz_one(B0);
    fmpz_one(B1);
    fmpz_one(order0);
    fmpz_one(order1);

    size_t t = s.size();
 
    for (int i = 1; i < t+1; i++) {
        if (i == 1) {
            fmpz_mul_ui(B1, B1, Primes[i]);
            fmpz_mul_ui(order0, order0, Primes[i]);
        }
        else if (i == 2) {
            fmpz_mul_ui(B0, B0, Primes[i]);
            fmpz_mul_ui(order1, order1, Primes[i]);
        }
        else {
            fmpz_mul_ui(B0, B0, Primes[i]);
            fmpz_mul_ui(B1, B1, Primes[i]);
        }
    }
    return { *B0, *B1, *order0, *order1 };
}

std::tuple<fmpz_t*, KummerLine, KummerPoint, KummerPoint> bin_terSIDH__hybrid::keygenA(flint_rand_t& prng, fmpz_t& A, fmpz_t& B, KummerPoint& xQA, KummerPoint& xPA, KummerPoint& xPQA, KummerPoint& xPB, KummerPoint& xQB, KummerLine& L0) {
    // generates secret and public keys for A

    // Generate the secret data
    fmpz_t skA;
    fmpz_randm(skA, prng, A);

    // compute the isogeny kernel
    KummerPoint xK = xQA.ladder_3_pt(xPA, xPQA, skA);

    // compute the isogeny
    KummerLineIsogenyComposite phiA = KummerLineIsogenyComposite(L0, xK, A);
    KummerLine EA = phiA.i_codomain;

    // generate the masking values
    fmpz_t mask, test, mask_inv;
    fmpz_set(mask, B);
    fmpz_gcd(test, mask, B);
    while (fmpz_is_one(test)) {
        fmpz_randm(test, prng, B + 1);
        fmpz_gcd(test, mask, B);
    }

    fmpz_invmod(mask_inv, mask, B);

    // scale the torsion images
    KummerPoint xR = phiA.call(xPB).multiply(mask);
    KummerPoint xS = phiA.call(xQB).multiply(mask_inv);

    return std::make_tuple(&skA, EA, xR, xS);
}

std::tuple<std::array<fmpz_t, 4>, KummerLine, KummerPoint, KummerPoint, KummerPoint> bin_terSIDH__hybrid::keygenB(std::vector<int> sk_choices, int t, KummerPoint& xPB, KummerPoint& xQB, KummerLine& L0, KummerPoint& xPA, KummerPoint& xQA, KummerPoint& xPQA) {
    // generates secret and public keys for B

    // generate the secret data
    std::vector<int> sk_calc;
    int idx;
    int size = int(sk_choices.size());
    for (int i = 0; i < t; ++i) {
        idx = std::rand() % size;
        sk_calc.push_back(sk_choices[idx]);
    }
    std::array<fmpz_t, 4> sk = compute_kernel_scalars_hybrid(sk_calc);

    // compute the isogeny kernels
    KummerPoint xK0 = xPB.multiply(sk[0]);
    KummerPoint xK1 = xQB.multiply(sk[1]);

    // compute the first isogeny
    KummerLineIsogenyComposite phiB0 = KummerLineIsogenyComposite(L0, xK0, sk[2]);
    xK1 = phiB0.call(xK1);

    // evaluate the first isogeny
    KummerPoint xPA0 = phiB0.call(xPA);
    KummerPoint xQA0 = phiB0.call(xQA);
    KummerPoint xPQA0 = phiB0.call(xPQA);

    // compute the second isogeny
    KummerLine EA0 = phiB0.i_codomain;
    KummerLineIsogenyComposite phiB1 = KummerLineIsogenyComposite(EA0, xK1, sk[3]);

    return std::make_tuple(sk, phiB1.i_codomain, phiB1.call(xPA0), phiB1.call(xQA0), phiB1.call(xPQA0));
}

fq_poly_t* bin_terSIDH__hybrid::sharedA(fmpz_t& skA, KummerLine& EB, KummerPoint& RA, KummerPoint& SA, KummerPoint& RSA, fmpz_t& A) {
    // compute shared j invariant from B's public key and A's private key

    // compute the isogeny kernel
    KummerPoint xK = SA.ladder_3_pt(RA, RSA, skA);

    // compute the isogeny
    KummerLineIsogenyComposite phiAdash = KummerLineIsogenyComposite(EB, xK, A);

    return phiAdash.i_codomain.curve_j_invariant();
}

fq_poly_t* bin_terSIDH__hybrid::sharedB(std::array<fmpz_t, 4> skB, KummerLine& EA, KummerPoint& RB, KummerPoint& SB) {
    // compute shared j invariant from A's public key and B's private key

    // compute the isogeny kernels
    KummerPoint xK0 = RB.multiply(skB[0]);
    KummerPoint xK1 = SB.multiply(skB[1]);

    // compute the first isogeny
    KummerLineIsogenyComposite phiBdash0 = KummerLineIsogenyComposite(EA, xK0, skB[2]);
    xK1 = phiBdash0.call(xK1);

    // compute the second isogeny
    KummerLine EAB0 = phiBdash0.i_codomain;
    KummerLineIsogenyComposite phiBdash1 = KummerLineIsogenyComposite(EAB0, xK1, skB[3]);

    return phiBdash1.i_codomain.curve_j_invariant();
}


void bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(int t, std::vector<int> sk_choices, int l) {
    // sets up and runs the bin and terSIDH hybrid protocol

    fmpz_t A, temp;
    fmpz_t* B;
    fmpz_set_ui(A, l);
    fmpz_mul_ui(A, A, 2);
    fmpz_set_ui(temp, 2);
    fmpz_pow_fmpz(A, temp, A);

    std::vector<slong> primes;
    for (int i = 1; i < t + 1; i++) {
        primes.push_back(Primes[i]);
    }

    B = ProductTree::prod(primes);

    fmpz_t p;
    fmpz_mul(p, A, *B);
    fmpz_print(p);

    std::cout << std::endl << "making primes" << std::endl << std::chrono::system_clock::now() << std::endl;
    int f = bin_terSIDH::make_prime(p, 1);
    std::cout << std::endl << "made primes" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "Size of field: ";
    fmpz_print(p);

    fq_ctx_t F;
    fmpz_mod_ctx_t FF;
    fmpz_mod_ctx_init(FF, p);
    fmpz_mod_poly_t x;
    fmpz_mod_poly_init(x, FF);
    fmpz_mod_poly_one(x, FF);
    fmpz_mod_poly_set_coeff_ui(x, 2, 1, FF);
    fq_ctx_init_modulus(F, x, FF, "i");
    fmpz_t order;
    fmpz_add_ui(order, p, 1);
    fmpz_pow_ui(order, order, 2);

    std::cout << std::endl << "making elliptic curve";
    EllipticCurve E0 = EllipticCurve(F, { 0, 6, 0, 1, 0 }, order, p);
    std::cout << std::endl << "made elliptic curve" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "Order of elliptic curve: ";
    fmpz_print(order);

    flint_rand_t prng;
    flint_randinit(prng);

    std::array<Point, 2> PAQA = TorsionData::torsion_basis(E0, A, p, prng);
    std::cout << std::endl << "generated torsion basis a" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::array<Point, 2> PBQB = TorsionData::torsion_basis(E0, *B, p, prng);
    std::cout << "generated torsion basis b" << std::endl << std::chrono::system_clock::now() << std::endl;

    fmpz_t A2;
    fmpz_divexact_ui(A, A, 2);
    Point test = Point::add_points(PAQA[0], PAQA[1]);

    // Ensures that 2*PA != (0, 0) and 2*QA != (0, 0) which causes problems with x only arithmetic
    if (fq_poly_is_zero(Point::scalar_multiplication(PAQA[0], A2).p_y, F)) {
        PAQA[0] = Point::add_points(PAQA[0], PAQA[1]);
    }
    else if (fq_poly_is_zero(Point::scalar_multiplication(test, A2).p_y, F)) {
        PAQA[1] = Point::add_points(PAQA[0], PAQA[1]);
    }

    test = Point::add_points(PAQA[0], PAQA[1]);

    if (fq_poly_is_zero(Point::scalar_multiplication(PAQA[0], A2).p_y, F)) {
        PAQA[0] = Point::add_points(PAQA[0], PAQA[1]);
    }
    else if (fq_poly_is_zero(Point::scalar_multiplication(test, A2).p_y, F)) {
        PAQA[1] = Point::add_points(PAQA[0], PAQA[1]);
    }

    Point PQA = Point::subtract(PAQA[0], PAQA[1]);
    Point PQB = Point::subtract(PBQB[0], PBQB[1]);

    KummerLine L0 = KummerLine(E0);

    KummerPoint xPA = KummerPoint(L0, *PAQA[0].x());
    KummerPoint xQA = KummerPoint(L0, *PAQA[1].x());
    KummerPoint xPQA = KummerPoint(L0, *PQA.x());

    KummerPoint xPB = KummerPoint(L0, *PBQB[0].x());
    KummerPoint xQB = KummerPoint(L0, *PBQB[1].x());
    KummerPoint xPQB = KummerPoint(L0, *PQB.x());

    fq_poly_t* ssA;
    fq_poly_t* ssB;

    std::cout << "generating a keys" << std::endl << std::chrono::system_clock::now() << std::endl;
    std::tuple<fmpz_t*, KummerLine, KummerPoint, KummerPoint> a_keys = keygenA(prng, A, *B, xQA, xPA, xPQA, xPB, xQB, L0);
    std::cout << "generated a keys" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generating b keys" << std::endl << std::chrono::system_clock::now() << std::endl;
    std::tuple<std::array<fmpz_t, 4>, KummerLine, KummerPoint, KummerPoint, KummerPoint> b_keys = keygenB(sk_choices, t, xPB, xQB, L0, xPA, xQA, xPQA);
    std::cout << "generated b keys" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generating a shared key" << std::endl << std::chrono::system_clock::now() << std::endl;
    ssA = sharedA(*std::get<0>(a_keys), std::get<1>(b_keys), std::get<2>(b_keys), std::get<3>(b_keys), std::get<4>(b_keys), A);
    std::cout << "generated a shared key" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generating b shared key" << std::endl << std::chrono::system_clock::now() << std::endl;
    ssB = sharedB(std::get<0>(b_keys), std::get<1>(a_keys), std::get<2>(a_keys), std::get<3>(a_keys));
    std::cout << "generated b shared key" << std::endl << std::chrono::system_clock::now() << std::endl;

    // checks shared key computed by A and B are equal
    bool answer = false;
    if (fq_poly_equal(*ssA, *ssB, F)) {
        answer = true;
    }

    std::cout << answer << std::endl;
}