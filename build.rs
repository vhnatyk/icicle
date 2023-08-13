use std::env;

fn main() {
    //TODO: check cargo features selected
    //TODO: can conflict/duplicate with make ?

    println!("cargo:rerun-if-env-changed=CXXFLAGS");
    println!("cargo:rerun-if-changed=./icicle");

    let arch_type = env::var("ARCH_TYPE").unwrap_or(String::from("sm_86"));
    let stream_type = env::var("DEFAULT_STREAM").unwrap_or(String::from("legacy"));

    let mut arch = String::from("-arch=");
    arch.push_str(&arch_type);
    let mut stream = String::from("-default-stream=");
    stream.push_str(&stream_type);

    let mut nvcc = cc::Build::new();

    println!("Compiling icicle library using arch: {}", &arch);

    if cfg!(feature = "g2") {
        nvcc.define("G2_DEFINED", None);
    }
    nvcc.cuda(true);
    nvcc.debug(false);
    nvcc.flag(&arch);
    nvcc.flag(&stream);
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
