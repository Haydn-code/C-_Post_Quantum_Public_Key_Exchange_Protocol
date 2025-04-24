#ifndef TESTING_H
#define TESTING_H

#include "bin-terSIDH--hybrid.h"

class testing {
	// class for testing implementations of bin and terSIDH protocols and their hybrid variants
public:
	static void core_testing();
	static void bin_terSIDH_tests();
	static void bin_terSIDH_hybrid_tests();
};

#endif