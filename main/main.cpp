#include "headers/main.h"

int main(int argc, char* argv[]) {
    // setup of security protocol and parameters with runtime arguments

    srand(42);

    int security = 128;
    int t = 0;
    bool ternary = true;
    bool hybrid = false;
    bool test = false;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--192") == 0) {
            security = 192;
            break;
        }
        if (std::strcmp(argv[i], "--256") == 0) {
            security = 256;
            break;
        }
    }

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--bin") == 0 or std::strcmp(argv[i], "--binary") == 0 ) {
            ternary = false;
            break;
        }
    }

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--hyb") == 0 or std::strcmp(argv[i], "--hybrid") == 0) {
            hybrid = true;
            break;
        }
    }

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--test") == 0) {
            test = true;
            break;
        }
    }

    int params_binSIDH[3] = { 134, 192, 256 };
    int params_terSIDH[3] = { 93, 128, 162 };
    int* params = params_binSIDH;
    std::vector<int> sk_choices;

    if (ternary) {
        params = params_terSIDH;
        sk_choices = { 0, 1, 2 };
    }
    else {
        sk_choices = { 1, 2 };
    }

    if (security == 128) {
        t = params[0];
    }
    else if (security == 192) {
        t = params[1];
    }
    else {
        t = params[2];
    }
    
    // either tests the implementation, or carries out the chosen protocol with the chosen parameters
    if (test) {
        testing::core_testing();
        return 0;
    }

    if (hybrid) {
        bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(t, sk_choices, security);
    }
    else {
        bin_terSIDH::bin_terSIDH_core(t, sk_choices);
    }
    return 0;
}
