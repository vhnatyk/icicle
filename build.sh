start_time=`date +%s`
nvcc -arch=sm_86 -maxrregcount=42 --threads=1280 --compiler-options=-pipe -o ingo_icicle ./icicle/curves/bls12_381/projective.cu ./icicle/curves/bls12_381/lde.cu ./icicle/curves/bls12_381/ve_mod_mult.cu && echo run time is $(expr `date +%s` - $start_time) s
