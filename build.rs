use std::env;

fn main() {
    //TODO: check cargo features selected
    //TODO: can conflict/duplicate with make ?

    println!("cargo:rerun-if-env-changed=CXXFLAGS");
    println!("cargo:rerun-if-changed=./icicle");

    let arch_type = env::var("ARCH_TYPE").unwrap_or(String::from("sm_86"));

    let mut arch = String::from("-arch=");
    arch.push_str(&arch_type);

    let mut nvcc = cc::Build::new();

    println!("Compiling icicle library using arch: {}", &arch);

    nvcc.cuda(true);
    nvcc.debug(false);
    nvcc.flag(&arch);
    // nvcc.flag("-maxrregcount=42");
    nvcc.flag("-w");
    nvcc.flag("--threads=1280");
    nvcc.flag("--compiler-options=-pipe");
    nvcc.files([    
        "./icicle/curves/bls12_381/projective.cu",
        "./icicle/curves/bls12_381/lde.cu",
        // "./icicle/curves/bls12_381/msm.cu",
        "./icicle/curves/bls12_381/ve_mod_mult.cu",
    ]);
    nvcc.compile("ingo_icicle"); //TODO: extension??
}
