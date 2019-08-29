# rvs-dft
Real-Valued Sequences DFT algorithms


* Computation of the N-point DFT of two N-point real-valued sequences
    - src/matlab/dft_split.m - Matlab algorithm
    - src/matlab/dft_split_tb.m - Matlab testbench
    - src/dft_split.c - C implementation


* Computation of the N/2-point DFT of a N-point real-valued sequence
    - src/matlab/dft_half.m - Matlab algorithm
    - src/matlab/dft_half_tb.m - Matlab testbench
    - src/dft_half.c - C implementation


* Computation of fast convolution of two N-point real-valued sequences with N/2-point DFT (type 1)
    - src/matlab/dft_fconv1.m - Matlab algorithm
    - src/matlab/dft_fconv1_tb.m - Matlab testbench
    - src/dft_fconv1.c - C implementation


* Computation of fast convolution of two N-point real-valued sequences with N/2-point DFT (type 2)
    - src/matlab/dft_fconv2.m - Matlab algorithm
    - src/matlab/dft_fconv2_tb.m - Matlab testbench
    - src/dft_fconv2.c - C implementation