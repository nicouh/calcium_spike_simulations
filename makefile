g++ main_ca_sims.cpp -O3 --std=c++17 -lstdc++fs -o caspikesims.out
echo "Arguments are <PosFdb RecoveryMode Runs dt MaxState Lambda rr_nrec rh_coef Nt g sp n Ka>"
echo "Example run:"
echo "./caspikesims.out exp normal 5 1e-05 10 0.005 1.23 00 10 2.50 1.0 3 00"


