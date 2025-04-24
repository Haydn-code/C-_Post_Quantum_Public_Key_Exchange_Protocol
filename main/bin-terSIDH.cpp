#include "headers/bin-terSIDH.h"

int bin_terSIDH::make_prime(fmpz_t& p, int f) {
    // given a value `p`, finds a cofactor `f` such that p*f -1 is prime

    fmpz_t temp;
    for (f; f < 1000; f++) {
        fmpz_mul_ui(temp, p, f);
        fmpz_sub_ui(temp, temp, 1);
        if (fmpz_is_prime(temp)) {
            fmpz_init_set(p, temp);
            fmpz_clear(temp);
            return f;
        }
    }
    throw std::invalid_argument("failed to generate primes");
}

std::array<fmpz_t, 4> bin_terSIDH::compute_kernel_scalars(std::vector<int> s, bool Alice) {
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
    std::array<int, 256>* P;

    if (Alice) {
        P = &Alice_Primes;
    }
    else {
        P = &Bob_Primes;
    }
    int p;
    for (int i = 0; i < 2*t; i++) {
        if (Alice && (i == 0)) {
            p = 4;
        }
        else {
            p = s[i];
        }
        if (s[i] == 2) {
            fmpz_mul_ui(B1, B1, p);
            fmpz_mul_ui(order0, order0, p);
        }
        else if (s[i] == 1) {
            fmpz_mul_ui(B0, B0, p);
            fmpz_mul_ui(order1, order1, p);
        }
        else {
            fmpz_mul_ui(B0, B0, p);
            fmpz_mul_ui(B1, B1, p);
        }
    }
    return { *B0, *B1, *order0, *order1 };
}

std::tuple<std::array<fmpz_t, 4>, KummerLine, std::array<KummerPoint, 2>> bin_terSIDH::keygen(std::array<KummerPoint, 2> A_points, std::array<KummerPoint, 2> B_points, std::vector<int> sk_choices, int t, KummerLine& L0, fmpz_t& modulus, flint_rand_t& prng, bool Alice) {
    // generates keys a public and private key pairing

    // Extract the Kummer points
    KummerPoint xP = A_points[0];
    KummerPoint xQ = A_points[1];
    KummerPoint xR = B_points[0];
    KummerPoint xS = B_points[1];

    //generate the secret data
    std::vector<int> sk_calc;
    int idx;
    int size = int(sk_choices.size());
    for (int i = 0; i < t; ++i) {
        idx = std::rand() % size;
        sk_calc.push_back(sk_choices[idx]);
    }
    std::array<fmpz_t, 4> sk = compute_kernel_scalars(sk_calc, Alice);

    //compute the isogeny kernels
    KummerPoint xK0 = xP.multiply(sk[0]);
    KummerPoint xK1 = xQ.multiply(sk[1]);

    //compute the first isogeny
    KummerLineIsogenyComposite phi0 = KummerLineIsogenyComposite(L0, xK0, sk[2]);
    xK1 = phi0.call(xK1);

    //evaluate action on aux data
    xR = phi0.call(xR);
    xS = phi0.call(xS);

    //compute the second isogeny from the codomain of phi0
    KummerLine E1 = phi0.i_codomain;
    KummerLineIsogenyComposite phi1 = KummerLineIsogenyComposite(E1, xK1, sk[3]);

    //evaluate the action on aux data
    xR = phi1.call(xR);
    xS = phi1.call(xS);

    fmpz_t mask, mask_inv, test;
    fmpz_set(mask, modulus);
    fmpz_gcd(test, mask, modulus);
    while (fmpz_is_one(test)) {
        fmpz_randm(mask, prng, modulus);
        fmpz_gcd(test, mask, modulus);
    }

    fmpz_invmod(mask_inv, mask, modulus);
    
    //scale the torsion images
    xR = xR.multiply(mask);
    xS = xS.multiply(mask_inv);
    std::array<KummerPoint, 2> xRxS = { xR, xS };

    //compute the public data
    E1 = phi1.i_codomain;

    return std::make_tuple(sk, E1, xRxS);
}

fq_poly_t* bin_terSIDH::shared(std::array<fmpz_t, 4> sk, KummerLine E, std::array<KummerPoint, 2> xRxS) {
    // computes the shared key from a public key and a private key

    // compute the isogeny kernels
    KummerPoint xK0 = xRxS[0].multiply(sk[0]);
    KummerPoint xK1 = xRxS[1].multiply(sk[1]);

    //compute the first isogeny
    KummerLineIsogenyComposite phi0(E, xK0, sk[2]);
    xK1 = phi0.call(xK1);

    //compute the second isogeny
    E = phi0.i_codomain;
    KummerLineIsogenyComposite phi1 = KummerLineIsogenyComposite(E, xK1, sk[3]);

    KummerLine EAB = phi1.i_codomain;

    return EAB.curve_j_invariant();
}

void bin_terSIDH::bin_terSIDH_core(int t, std::vector<int> sk_choices) {
    // sets up bin_terSIDH protocols and runs key exchange of selected protocol.
    
    // setup
    fmpz_t A, B;
    fmpz_init_set_si(A, 1);
    fmpz_init_set_si(B, 1);

    for (int i = 0; i < 2 * t; i++) {
        if (i % 2 == 0) {
            fmpz_mul_ui(A, A, Primes[i]);
        }
        else {
            fmpz_mul_ui(B, B, Primes[i]);
        }
    }
    fmpz_mul_si(A, A, 2);

    std::cout << "Value of A: ";
    fmpz_print(A);
    std::cout << std::endl << "Value of B: ";
    fmpz_print(B);
    fmpz_t p;
    fmpz_mul(p, A, B);
    
    std::cout << std::endl << "making primes" << std::endl << std::chrono::system_clock::now() << std::endl;
    int f = make_prime(p, 1);
    std::cout << "made primes" << std::endl << std::chrono::system_clock::now() << std::endl;

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
    
    std::cout << std::endl << "making elliptic curve" << std::endl << std::chrono::system_clock::now();
    EllipticCurve E0 = EllipticCurve(F, { 0, 6, 0, 1, 0 }, order, p);
    std::cout << std::endl << "made elliptic curve" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "Order of elliptic curve: ";
    fmpz_print(order);

    flint_rand_t prng;
    flint_randinit(prng);

    std::cout << std::endl << "generating torsion basis a" << std::endl << std::chrono::system_clock::now() << std::endl;
    std::array<Point, 2> PAQA = TorsionData::torsion_basis(E0, A, p, prng);
    std::cout << "generated torsion basis a" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generating torsion basis b" << std::endl << std::chrono::system_clock::now() << std::endl;
    std::array<Point, 2> PBQB = TorsionData::torsion_basis(E0, B, p, prng);
    std::cout << "generated torsion basis b" << std::endl << std::chrono::system_clock::now() << std::endl;
    
    
    fmpz_t A2;
    fmpz_divexact_ui(A, A, 2);

    // Ensures that 2*PA != (0, 0) and 2*QA != (0, 0) which causes problems with x only arithmetic
    if (fq_poly_is_zero(Point::scalar_multiplication(PAQA[0], A2).p_y, F)) {
        PAQA[0] = Point::add_points(PAQA[0], PAQA[1]);
    }
    else if (fq_poly_is_zero(Point::scalar_multiplication(PAQA[1], A2).p_y, F)) {
        PAQA[1] = Point::add_points(PAQA[0], PAQA[1]);
    }

    KummerLine L0 = KummerLine(E0);

    KummerPoint xPA = KummerPoint(L0, *PAQA[0].x());
    KummerPoint xQA = KummerPoint(L0, *PAQA[1].x());

    KummerPoint xPB = KummerPoint(L0, *PBQB[0].x());
    KummerPoint xQB = KummerPoint(L0, *PBQB[1].x());

    std::cout << "generating a keys" << std::endl << std::chrono::system_clock::now() << std::endl;
    std::tuple<std::array<fmpz_t, 4>, KummerLine, std::array<KummerPoint, 2>> a_keys = keygen({xPA, xQA}, {xPB, xQB}, sk_choices, t, L0, B, prng, true);
    std::cout << "generated a keys" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generated b keys" << std::endl << std::chrono::system_clock::now() << std::endl;
    std::tuple<std::array<fmpz_t, 4>, KummerLine, std::array<KummerPoint, 2>> b_keys = keygen({ xPB, xQB }, { xPA, xQA }, sk_choices, t, L0, A, prng, false);
    std::cout << "generated b keys" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generating a shared key" << std::endl << std::chrono::system_clock::now() << std::endl;
    fq_poly_t* ssA = shared(std::get<0>(a_keys), std::get<1>(a_keys), std::get<2>(a_keys));
    std::cout << "generated a shared key" << std::endl << std::chrono::system_clock::now() << std::endl;

    std::cout << "generating b shared key" << std::endl << std::chrono::system_clock::now() << std::endl;
    fq_poly_t* ssB = shared(std::get<0>(b_keys), std::get<1>(b_keys), std::get<2>(b_keys));
    std::cout << "generated b shared key" << std::endl << std::chrono::system_clock::now() << std::endl;

    // check shared key computed by A and B are equal
    bool test = false;
    if (fq_poly_equal(*ssA, *ssB, F)) {
        test = true;
    }

    std::cout << test << std::endl;
}