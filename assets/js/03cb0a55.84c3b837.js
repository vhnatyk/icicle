"use strict";(self.webpackChunkdocusaurus=self.webpackChunkdocusaurus||[]).push([[912],{5680:(e,n,t)=>{t.d(n,{xA:()=>u,yg:()=>g});var i=t(6540);function a(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function o(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);n&&(i=i.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,i)}return t}function r(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?o(Object(t),!0).forEach((function(n){a(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):o(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,i,a=function(e,n){if(null==e)return{};var t,i,a={},o=Object.keys(e);for(i=0;i<o.length;i++)t=o[i],n.indexOf(t)>=0||(a[t]=e[t]);return a}(e,n);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(i=0;i<o.length;i++)t=o[i],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(a[t]=e[t])}return a}var s=i.createContext({}),p=function(e){var n=i.useContext(s),t=n;return e&&(t="function"==typeof e?e(n):r(r({},n),e)),t},u=function(e){var n=p(e.components);return i.createElement(s.Provider,{value:n},e.children)},d="mdxType",c={inlineCode:"code",wrapper:function(e){var n=e.children;return i.createElement(i.Fragment,{},n)}},m=i.forwardRef((function(e,n){var t=e.components,a=e.mdxType,o=e.originalType,s=e.parentName,u=l(e,["components","mdxType","originalType","parentName"]),d=p(t),m=a,g=d["".concat(s,".").concat(m)]||d[m]||c[m]||o;return t?i.createElement(g,r(r({ref:n},u),{},{components:t})):i.createElement(g,r({ref:n},u))}));function g(e,n){var t=arguments,a=n&&n.mdxType;if("string"==typeof e||a){var o=t.length,r=new Array(o);r[0]=m;var l={};for(var s in n)hasOwnProperty.call(n,s)&&(l[s]=n[s]);l.originalType=e,l[d]="string"==typeof e?e:a,r[1]=l;for(var p=2;p<o;p++)r[p]=t[p];return i.createElement.apply(null,r)}return i.createElement.apply(null,t)}m.displayName="MDXCreateElement"},8542:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>s,contentTitle:()=>r,default:()=>c,frontMatter:()=>o,metadata:()=>l,toc:()=>p});var i=t(8168),a=(t(6540),t(5680));t(1873);const o={},r="Poseidon",l={unversionedId:"icicle/primitives/poseidon",id:"icicle/primitives/poseidon",title:"Poseidon",description:"Poseidon is a popular hash in the ZK ecosystem primarily because its optimized to work over large prime fields, a common setting for ZK proofs, thereby minimizing the number of multiplicative operations required.",source:"@site/docs/icicle/primitives/poseidon.md",sourceDirName:"icicle/primitives",slug:"/icicle/primitives/poseidon",permalink:"/icicle/primitives/poseidon",editUrl:"https://github.com/ingonyama-zk/icicle/tree/main/docs/icicle/primitives/poseidon.md",tags:[],version:"current",lastUpdatedBy:"VitaliiH",lastUpdatedAt:1715232285,formattedLastUpdatedAt:"5/9/2024",frontMatter:{},sidebar:"GettingStartedSidebar",previous:{title:"Keccak",permalink:"/icicle/primitives/keccak"},next:{title:"Polynomial API Overview",permalink:"/icicle/polynomials/overview"}},s={},p=[{value:"Initialization",id:"initialization",level:2},{value:"Applying full and partial rounds",id:"applying-full-and-partial-rounds",level:2},{value:"Full rounds",id:"full-rounds",level:3},{value:"Partial Rounds",id:"partial-rounds",level:3},{value:"Using Poseidon",id:"using-poseidon",level:2},{value:"Supported Bindings",id:"supported-bindings",level:3},{value:"Constants",id:"constants",level:3},{value:"Rust API",id:"rust-api",level:3},{value:"The Tree Builder",id:"the-tree-builder",level:2},{value:"Benchmarks",id:"benchmarks",level:3}],u={toc:p},d="wrapper";function c(e){let{components:n,...o}=e;return(0,a.yg)(d,(0,i.A)({},u,o,{components:n,mdxType:"MDXLayout"}),(0,a.yg)("h1",{id:"poseidon"},"Poseidon"),(0,a.yg)("p",null,(0,a.yg)("a",{parentName:"p",href:"https://eprint.iacr.org/2019/458.pdf"},"Poseidon")," is a popular hash in the ZK ecosystem primarily because its optimized to work over large prime fields, a common setting for ZK proofs, thereby minimizing the number of multiplicative operations required."),(0,a.yg)("p",null,"Poseidon has also been specifically designed to be efficient when implemented within ZK circuits, Poseidon uses far less constraints compared to other hash functions like Keccak or SHA-256 in the context of ZK circuits."),(0,a.yg)("p",null,"Poseidon has been used in many popular ZK protocols such as Filecoin and ",(0,a.yg)("a",{parentName:"p",href:"https://drive.google.com/file/d/1bZZvKMQHaZGA4L9eZhupQLyGINkkFG_b/view?usp=drive_open"},"Plonk"),"."),(0,a.yg)("p",null,"Our implementation of Poseidon is implemented in accordance with the optimized ",(0,a.yg)("a",{parentName:"p",href:"https://spec.filecoin.io/algorithms/crypto/poseidon/"},"Filecoin version"),"."),(0,a.yg)("p",null,"Lets understand how Poseidon works."),(0,a.yg)("h2",{id:"initialization"},"Initialization"),(0,a.yg)("p",null,"Poseidon starts with the initialization of its internal state, which is composed of the input elements and some pre-generated constants. An initial round constant is added to each element of the internal state. Adding the round constants ensures the state is properly mixed from the beginning."),(0,a.yg)("p",null,"This is done to prevent collisions and to prevent certain cryptographic attacks by ensuring that the internal state is sufficiently mixed and unpredictable."),(0,a.yg)("p",null,(0,a.yg)("img",{alt:"Alt text",src:t(9569).A,width:"872",height:"294"})),(0,a.yg)("h2",{id:"applying-full-and-partial-rounds"},"Applying full and partial rounds"),(0,a.yg)("p",null,'To generate a secure hash output, the algorithm goes through a series of "full rounds" and "partial rounds" as well as transformations between these sets of rounds in the following order:'),(0,a.yg)("p",null,(0,a.yg)("inlineCode",{parentName:"p"},"First full rounds -> apply S-box and Round constants -> partial rounds -> Last full rounds -> Apply S-box")),(0,a.yg)("h3",{id:"full-rounds"},"Full rounds"),(0,a.yg)("p",null,(0,a.yg)("img",{alt:"Alt text",src:t(2811).A,width:"864",height:"552"})),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"Uniform Application of S-box:")," In full rounds, the S-box (a non-linear transformation) is applied uniformly to every element of the hash function's internal state. This ensures a high degree of mixing and diffusion, contributing to the hash function's security. The functions S-box involves raising each element of the state to a certain power denoted by ",(0,a.yg)("inlineCode",{parentName:"p"},"\u03b1")," a member of the finite field defined by the prime ",(0,a.yg)("inlineCode",{parentName:"p"},"p"),"; ",(0,a.yg)("inlineCode",{parentName:"p"},"\u03b1")," can be different depending on the the implementation and user configuration."),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"Linear Transformation:")," After applying the S-box, a linear transformation is performed on the state. This involves multiplying the state by a MDS (Maximum Distance Separable) Matrix. which further diffuses the transformations applied by the S-box across the entire state."),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"Addition of Round Constants:")," Each element of the state is then modified by adding a unique round constant. These constants are different for each round and are precomputed as part of the hash function's initialization. The addition of round constants ensures that even minor changes to the input produce significant differences in the output."),(0,a.yg)("h3",{id:"partial-rounds"},"Partial Rounds"),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"Selective Application of S-Box:")," Partial rounds apply the S-box transformation to only one element of the internal state per round, rather than to all elements. This selective application significantly reduces the computational complexity of the hash function without compromising its security. The choice of which element to apply the S-box to can follow a specific pattern or be fixed, depending on the design of the hash function."),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"Linear Transformation and Round Constants:")," A linear transformation is performed and round constants are added. The linear transformation in partial rounds can be designed to be less computationally intensive (this is done by using a sparse matrix) than in full rounds, further optimizing the function's efficiency."),(0,a.yg)("p",null,"The user of Poseidon can often choose how many partial or full rounds he wishes to apply; more full rounds will increase security but degrade performance. The choice and balance is highly dependent on the use case."),(0,a.yg)("p",null,(0,a.yg)("img",{alt:"Alt text",src:t(7408).A,width:"866",height:"560"})),(0,a.yg)("h2",{id:"using-poseidon"},"Using Poseidon"),(0,a.yg)("p",null,"ICICLE Poseidon is implemented for GPU and parallelization is performed for each element of the state rather than for each state.\nWhat that means is we calculate multiple hash-sums over multiple pre-images in parallel, rather than going block by block over the input vector."),(0,a.yg)("p",null,"So for Poseidon of arity 2 and input of size 1024 * 2, we would expect 1024 elements of output. Which means each block would be of size 2 and that would result in 1024 Poseidon hashes being performed."),(0,a.yg)("h3",{id:"supported-bindings"},"Supported Bindings"),(0,a.yg)("p",null,(0,a.yg)("a",{parentName:"p",href:"https://github.com/ingonyama-zk/icicle/tree/main/wrappers/rust/icicle-core/src/poseidon"},(0,a.yg)("inlineCode",{parentName:"a"},"Rust"))),(0,a.yg)("h3",{id:"constants"},"Constants"),(0,a.yg)("p",null,"Poseidon is extremely customizable and using different constants will produce different hashes, security levels and performance results."),(0,a.yg)("p",null,"We support pre-calculated and optimized constants for each of the ",(0,a.yg)("a",{parentName:"p",href:"#supported-curves"},"supported curves"),".The constants can be found ",(0,a.yg)("a",{parentName:"p",href:"https://github.com/ingonyama-zk/icicle/tree/main/icicle/include/poseidon/constants"},"here")," and are labeled clearly per curve ",(0,a.yg)("inlineCode",{parentName:"p"},"<curve_name>_poseidon.h"),"."),(0,a.yg)("p",null,"If you wish to generate your own constants you can use our python script which can be found ",(0,a.yg)("a",{parentName:"p",href:"https://github.com/ingonyama-zk/icicle/tree/main/icicle/include/poseidon/constants/generate_parameters.py"},"here"),"."),(0,a.yg)("p",null,"Prerequisites:"),(0,a.yg)("ul",null,(0,a.yg)("li",{parentName:"ul"},"Install python 3"),(0,a.yg)("li",{parentName:"ul"},(0,a.yg)("inlineCode",{parentName:"li"},"pip install poseidon-hash")),(0,a.yg)("li",{parentName:"ul"},(0,a.yg)("inlineCode",{parentName:"li"},"pip install galois==0.3.7")),(0,a.yg)("li",{parentName:"ul"},(0,a.yg)("inlineCode",{parentName:"li"},"pip install numpy"))),(0,a.yg)("p",null,"You will then need to modify the following values before running the script."),(0,a.yg)("pre",null,(0,a.yg)("code",{parentName:"pre",className:"language-python"},"# Modify these\narity = 11 # we support arity 2, 4, 8 and 11.\np = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001 # bls12-381\n# p = 0x12ab655e9a2ca55660b44d1e5c37b00159aa76fed00000010a11800000000001 # bls12-377\n# p = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001 # bn254\n# p = 0x1ae3a4617c510eac63b05c06ca1493b1a22d9f300f5138f1ef3622fba094800170b5d44300000008508c00000000001 # bw6-761\nprime_bit_len = 255\nfield_bytes = 32\n\n...\n\n# primitive_element = None\nprimitive_element = 7 # bls12-381\n# primitive_element = 22 # bls12-377\n# primitive_element = 5 # bn254\n# primitive_element = 15 # bw6-761\n")),(0,a.yg)("p",null,"We only support ",(0,a.yg)("inlineCode",{parentName:"p"},"alpha = 5")," so if you want to use another alpha for S-box please reach out on discord or open a github issue."),(0,a.yg)("h3",{id:"rust-api"},"Rust API"),(0,a.yg)("p",null,"This is the most basic way to use the Poseidon API."),(0,a.yg)("pre",null,(0,a.yg)("code",{parentName:"pre",className:"language-rust"},"let test_size = 1 << 10;\nlet arity = 2u32;\nlet ctx = get_default_device_context();\nlet constants = load_optimized_poseidon_constants::<F>(arity, &ctx).unwrap();\nlet config = PoseidonConfig::default();\n\nlet inputs = vec![F::one(); test_size * arity as usize];\nlet outputs = vec![F::zero(); test_size];\nlet mut input_slice = HostOrDeviceSlice::on_host(inputs);\nlet mut output_slice = HostOrDeviceSlice::on_host(outputs);\n\nposeidon_hash_many::<F>(\n    &mut input_slice,\n    &mut output_slice,\n    test_size as u32,\n    arity as u32,\n    &constants,\n    &config,\n)\n.unwrap();\n")),(0,a.yg)("p",null,"The ",(0,a.yg)("inlineCode",{parentName:"p"},"PoseidonConfig::default()")," can be modified, by default the inputs and outputs are set to be on ",(0,a.yg)("inlineCode",{parentName:"p"},"Host")," for example."),(0,a.yg)("pre",null,(0,a.yg)("code",{parentName:"pre",className:"language-rust"},"impl<'a> Default for PoseidonConfig<'a> {\n    fn default() -> Self {\n        let ctx = get_default_device_context();\n        Self {\n            ctx,\n            are_inputs_on_device: false,\n            are_outputs_on_device: false,\n            input_is_a_state: false,\n            aligned: false,\n            loop_state: false,\n            is_async: false,\n        }\n    }\n}\n")),(0,a.yg)("p",null,"In the example above ",(0,a.yg)("inlineCode",{parentName:"p"},"load_optimized_poseidon_constants::<F>(arity, &ctx).unwrap();")," is used which will load the correct constants based on arity and curve. Its possible to ",(0,a.yg)("a",{parentName:"p",href:"#constants"},"generate")," your own constants and load them."),(0,a.yg)("pre",null,(0,a.yg)("code",{parentName:"pre",className:"language-rust"},'let ctx = get_default_device_context();\n    let cargo_manifest_dir = env!("CARGO_MANIFEST_DIR");\n    let constants_file = PathBuf::from(cargo_manifest_dir)\n        .join("tests")\n        .join(format!("{}_constants.bin", field_prefix));\n    let mut constants_buf = vec![];\n    File::open(constants_file)\n        .unwrap()\n        .read_to_end(&mut constants_buf)\n        .unwrap();\n\n    let mut custom_constants = vec![];\n    for chunk in constants_buf.chunks(field_bytes) {\n        custom_constants.push(F::from_bytes_le(chunk));\n    }\n\n    let custom_constants = create_optimized_poseidon_constants::<F>(\n        arity as u32,\n        &ctx,\n        full_rounds_half,\n        partial_rounds,\n        &mut custom_constants,\n    )\n    .unwrap();\n')),(0,a.yg)("h2",{id:"the-tree-builder"},"The Tree Builder"),(0,a.yg)("p",null,"The tree builder allows you to build Merkle trees using Poseidon."),(0,a.yg)("p",null,"You can define both the tree's ",(0,a.yg)("inlineCode",{parentName:"p"},"height")," and its ",(0,a.yg)("inlineCode",{parentName:"p"},"arity"),". The tree ",(0,a.yg)("inlineCode",{parentName:"p"},"height")," determines the number of layers in the tree, including the root and the leaf layer. The ",(0,a.yg)("inlineCode",{parentName:"p"},"arity")," determines how many children each internal node can have."),(0,a.yg)("pre",null,(0,a.yg)("code",{parentName:"pre",className:"language-rust"},'let height = 20;\nlet arity = 2;\nlet leaves = vec![F::one(); 1 << (height - 1)];\nlet mut digests = vec![F::zero(); merkle_tree_digests_len(height, arity)];\n\nlet mut leaves_slice = HostOrDeviceSlice::on_host(leaves);\n\nlet ctx = get_default_device_context();\nlet constants = load_optimized_poseidon_constants::<F>(arity, &ctx).unwrap()\n\nlet mut config = TreeBuilderConfig::default();\nconfig.keep_rows = 1;\nbuild_poseidon_merkle_tree::<F>(&mut leaves_slice, &mut digests, height, arity, &constants, &config).unwrap();\n\nprintln!("Root: {:?}", digests[0..1][0]);\n')),(0,a.yg)("p",null,"Similar to Poseidon, you can also configure the Tree Builder ",(0,a.yg)("inlineCode",{parentName:"p"},"TreeBuilderConfig::default()")),(0,a.yg)("ul",null,(0,a.yg)("li",{parentName:"ul"},(0,a.yg)("inlineCode",{parentName:"li"},"keep_rows"),": The number of rows which will be written to output, 0 will write all rows."),(0,a.yg)("li",{parentName:"ul"},(0,a.yg)("inlineCode",{parentName:"li"},"are_inputs_on_device"),": Have the inputs been loaded to device memory ?"),(0,a.yg)("li",{parentName:"ul"},(0,a.yg)("inlineCode",{parentName:"li"},"is_async"),": Should the TreeBuilder run asynchronously? ",(0,a.yg)("inlineCode",{parentName:"li"},"False")," will block the current CPU thread. ",(0,a.yg)("inlineCode",{parentName:"li"},"True")," will require you call ",(0,a.yg)("inlineCode",{parentName:"li"},"cudaStreamSynchronize")," or ",(0,a.yg)("inlineCode",{parentName:"li"},"cudaDeviceSynchronize")," to retrieve the result.")),(0,a.yg)("h3",{id:"benchmarks"},"Benchmarks"),(0,a.yg)("p",null,"We ran the Poseidon tree builder on:"),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"CPU"),": 12th Gen Intel(R) Core(TM) i9-12900K/"),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"GPU"),": RTX 3090 Ti"),(0,a.yg)("p",null,(0,a.yg)("strong",{parentName:"p"},"Tree height"),": 30 (2^29 elements)"),(0,a.yg)("p",null,"The benchmarks include copying data from and to the device."),(0,a.yg)("table",null,(0,a.yg)("thead",{parentName:"table"},(0,a.yg)("tr",{parentName:"thead"},(0,a.yg)("th",{parentName:"tr",align:null},"Rows to keep parameter"),(0,a.yg)("th",{parentName:"tr",align:null},"Run time, Icicle"),(0,a.yg)("th",{parentName:"tr",align:null},"Supranational PC2"))),(0,a.yg)("tbody",{parentName:"table"},(0,a.yg)("tr",{parentName:"tbody"},(0,a.yg)("td",{parentName:"tr",align:null},"10"),(0,a.yg)("td",{parentName:"tr",align:null},"9.4 seconds"),(0,a.yg)("td",{parentName:"tr",align:null},"13.6 seconds")),(0,a.yg)("tr",{parentName:"tbody"},(0,a.yg)("td",{parentName:"tr",align:null},"20"),(0,a.yg)("td",{parentName:"tr",align:null},"9.5 seconds"),(0,a.yg)("td",{parentName:"tr",align:null},"13.6 seconds")),(0,a.yg)("tr",{parentName:"tbody"},(0,a.yg)("td",{parentName:"tr",align:null},"29"),(0,a.yg)("td",{parentName:"tr",align:null},"13.7 seconds"),(0,a.yg)("td",{parentName:"tr",align:null},"13.6 seconds")))))}c.isMDXComponent=!0},2811:(e,n,t)=>{t.d(n,{A:()=>i});const i=t.p+"assets/images/image-1-5902a94e1680e802186934fdf8ff205e.png"},7408:(e,n,t)=>{t.d(n,{A:()=>i});const i=t.p+"assets/images/image-2-943ff9b12b39ca32378f75f981ebfb7b.png"},9569:(e,n,t)=>{t.d(n,{A:()=>i});const i=t.p+"assets/images/image-bff569244d897bcfcd16d979cb29fb9c.png"},1873:(e,n,t)=>{t(6540)}}]);