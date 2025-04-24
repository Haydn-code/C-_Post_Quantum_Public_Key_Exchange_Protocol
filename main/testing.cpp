#include "headers/testing.h"

void testing::core_testing() {
	// runs tests for binSIDH and terSIDH protocols and their hybrid variants

	bin_terSIDH_tests();
	bin_terSIDH_hybrid_tests();
}

void testing::bin_terSIDH_tests(){
	//tests for binSIDH and terSIDH protocols

	// tests for binSIDH protocol
	bin_terSIDH::bin_terSIDH_core(134, { 1, 2 });
	bin_terSIDH::bin_terSIDH_core(192, { 1, 2 });
	bin_terSIDH::bin_terSIDH_core(256, { 1, 2 });

	//tests for terSIDH protocol
	bin_terSIDH::bin_terSIDH_core(93, { 0, 1, 2 });
	bin_terSIDH::bin_terSIDH_core(128, { 0, 1, 2 });
	bin_terSIDH::bin_terSIDH_core(162, { 0, 1, 2 });
}

void testing::bin_terSIDH_hybrid_tests(){
	//tests for binSIDH and terSIDH hybrid protocols

	// tests for binSIDH hybrid protocol
	bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(134, { 1, 2 }, 128);
	bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(192, { 1, 2 }, 192);
	bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(256, { 1, 2 }, 256);

	//tests for terSIDH protocol protocol
	bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(93, { 0, 1, 2 }, 128);
	bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(128, { 0, 1, 2 }, 192);
	bin_terSIDH__hybrid::bin_terSIDH__hybrid_core(162, { 0, 1, 2 }, 256);
}
