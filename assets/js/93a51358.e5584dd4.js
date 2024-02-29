"use strict";(self.webpackChunkdocusaurus=self.webpackChunkdocusaurus||[]).push([[966],{5680:(e,n,i)=>{i.d(n,{xA:()=>u,yg:()=>m});var t=i(6540);function r(e,n,i){return n in e?Object.defineProperty(e,n,{value:i,enumerable:!0,configurable:!0,writable:!0}):e[n]=i,e}function a(e,n){var i=Object.keys(e);if(Object.getOwnPropertySymbols){var t=Object.getOwnPropertySymbols(e);n&&(t=t.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),i.push.apply(i,t)}return i}function c(e){for(var n=1;n<arguments.length;n++){var i=null!=arguments[n]?arguments[n]:{};n%2?a(Object(i),!0).forEach((function(n){r(e,n,i[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(i)):a(Object(i)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(i,n))}))}return e}function o(e,n){if(null==e)return{};var i,t,r=function(e,n){if(null==e)return{};var i,t,r={},a=Object.keys(e);for(t=0;t<a.length;t++)i=a[t],n.indexOf(i)>=0||(r[i]=e[i]);return r}(e,n);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(t=0;t<a.length;t++)i=a[t],n.indexOf(i)>=0||Object.prototype.propertyIsEnumerable.call(e,i)&&(r[i]=e[i])}return r}var l=t.createContext({}),s=function(e){var n=t.useContext(l),i=n;return e&&(i="function"==typeof e?e(n):c(c({},n),e)),i},u=function(e){var n=s(e.components);return t.createElement(l.Provider,{value:n},e.children)},p="mdxType",d={inlineCode:"code",wrapper:function(e){var n=e.children;return t.createElement(t.Fragment,{},n)}},g=t.forwardRef((function(e,n){var i=e.components,r=e.mdxType,a=e.originalType,l=e.parentName,u=o(e,["components","mdxType","originalType","parentName"]),p=s(i),g=r,m=p["".concat(l,".").concat(g)]||p[g]||d[g]||a;return i?t.createElement(m,c(c({ref:n},u),{},{components:i})):t.createElement(m,c({ref:n},u))}));function m(e,n){var i=arguments,r=n&&n.mdxType;if("string"==typeof e||r){var a=i.length,c=new Array(a);c[0]=g;var o={};for(var l in n)hasOwnProperty.call(n,l)&&(o[l]=n[l]);o.originalType=e,o[p]="string"==typeof e?e:r,c[1]=o;for(var s=2;s<a;s++)c[s]=i[s];return t.createElement.apply(null,c)}return t.createElement.apply(null,i)}g.displayName="MDXCreateElement"},3703:(e,n,i)=>{i.r(n),i.d(n,{assets:()=>l,contentTitle:()=>c,default:()=>d,frontMatter:()=>a,metadata:()=>o,toc:()=>s});var t=i(8168),r=(i(6540),i(5680));i(1873);const a={},c="Rust bindings",o={unversionedId:"icicle/rust-bindings",id:"icicle/rust-bindings",title:"Rust bindings",description:"Rust bindings allow you to use ICICLE as a rust library.",source:"@site/docs/icicle/rust-bindings.md",sourceDirName:"icicle",slug:"/icicle/rust-bindings",permalink:"/icicle/icicle/rust-bindings",editUrl:"https://github.com/ingonyama-zk/developer-docs/tree/main/docs/icicle/rust-bindings.md",tags:[],version:"current",lastUpdatedBy:"Jeremy Felder",lastUpdatedAt:1709188030,formattedLastUpdatedAt:"2/29/2024",frontMatter:{},sidebar:"GettingStartedSidebar",previous:{title:"Golang bindings",permalink:"/icicle/icicle/golang-bindings"},next:{title:"Multi GPU APIs",permalink:"/icicle/icicle/rust-bindings/multi-gpu"}},l={},s=[{value:"Using ICICLE Rust bindings in your project",id:"using-icicle-rust-bindings-in-your-project",level:2}],u={toc:s},p="wrapper";function d(e){let{components:n,...i}=e;return(0,r.yg)(p,(0,t.A)({},u,i,{components:n,mdxType:"MDXLayout"}),(0,r.yg)("h1",{id:"rust-bindings"},"Rust bindings"),(0,r.yg)("p",null,"Rust bindings allow you to use ICICLE as a rust library."),(0,r.yg)("p",null,(0,r.yg)("inlineCode",{parentName:"p"},"icicle-core")," defines all interfaces, macros and common methods."),(0,r.yg)("p",null,(0,r.yg)("inlineCode",{parentName:"p"},"icicle-cuda-runtime")," defines DeviceContext which can be used to manage a specific GPU as well as wrapping common CUDA methods."),(0,r.yg)("p",null,(0,r.yg)("inlineCode",{parentName:"p"},"icicle-curves")," implements all interfaces and macros from icicle-core for each curve. For example icicle-bn254 implements curve bn254. Each curve has its own build script which will build the CUDA libraries for that curve as part of the rust-toolchain build."),(0,r.yg)("h2",{id:"using-icicle-rust-bindings-in-your-project"},"Using ICICLE Rust bindings in your project"),(0,r.yg)("p",null,"Simply add the following to your ",(0,r.yg)("inlineCode",{parentName:"p"},"Cargo.toml"),"."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre"},'# GPU Icicle integration\nicicle-cuda-runtime = { git = "https://github.com/ingonyama-zk/icicle.git" }\nicicle-core = { git = "https://github.com/ingonyama-zk/icicle.git" }\nicicle-bn254 = { git = "https://github.com/ingonyama-zk/icicle.git" }\n')),(0,r.yg)("p",null,(0,r.yg)("inlineCode",{parentName:"p"},"icicle-bn254")," being the curve you wish to use and ",(0,r.yg)("inlineCode",{parentName:"p"},"icicle-core")," and ",(0,r.yg)("inlineCode",{parentName:"p"},"icicle-cuda-runtime")," contain ICICLE utilities and CUDA wrappers."),(0,r.yg)("p",null,"If you wish to point to a specific ICICLE branch add ",(0,r.yg)("inlineCode",{parentName:"p"},'branch = "<name_of_branch>"')," or ",(0,r.yg)("inlineCode",{parentName:"p"},'tag = "<name_of_tag>"')," to the ICICLE dependency. For a specific commit add ",(0,r.yg)("inlineCode",{parentName:"p"},'rev = "<commit_id>"'),"."),(0,r.yg)("p",null,"When you build your project ICICLE will be built as part of the build command."),(0,r.yg)("h1",{id:"how-do-the-rust-bindings-work"},"How do the rust bindings work?"),(0,r.yg)("p",null,"The rust bindings are just rust wrappers for ICICLE Core static libraries which can be compiled. We integrate the compilation of the static libraries into rusts toolchain to make usage seamless and easy. This is achieved by ",(0,r.yg)("a",{parentName:"p",href:"https://github.com/ingonyama-zk/icicle/blob/main/wrappers/rust/icicle-curves/icicle-bn254/build.rs"},"extending rusts build command"),"."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-rust"},'use cmake::Config;\nuse std::env::var;\n\nfn main() {\n    println!("cargo:rerun-if-env-changed=CXXFLAGS");\n    println!("cargo:rerun-if-changed=../../../../icicle");\n\n    let cargo_dir = var("CARGO_MANIFEST_DIR").unwrap();\n    let profile = var("PROFILE").unwrap();\n\n    let out_dir = Config::new("../../../../icicle")\n                .define("BUILD_TESTS", "OFF") //TODO: feature\n                .define("CURVE", "bn254")\n                .define("CMAKE_BUILD_TYPE", "Release")\n                .build_target("icicle")\n                .build();\n\n    println!("cargo:rustc-link-search={}/build", out_dir.display());\n\n    println!("cargo:rustc-link-lib=ingo_bn254");\n    println!("cargo:rustc-link-lib=stdc++");\n    // println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");\n    println!("cargo:rustc-link-lib=cudart");\n}\n')))}d.isMDXComponent=!0},1873:(e,n,i)=>{i(6540)}}]);