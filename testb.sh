

flag='export RUSTFLAGS=""'
if [[ $2 = "f" ]]; then   
flag="export RUSTFLAGS=-Awarnings"
fi

watch=""
if [[ $3 = "w" ]]; then   
watch="watch -n 10"
fi

CMD="$watch cargo +nightly test --release --package icicle-utils --lib -- test_bls12_381::tests_bls12_381::$1 --exact --nocapture"

echo cmd is: $CMD
$flag; $CMD;
