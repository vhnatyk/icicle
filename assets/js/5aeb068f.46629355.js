"use strict";(self.webpackChunkdocusaurus=self.webpackChunkdocusaurus||[]).push([[544],{5680:(e,n,t)=>{t.d(n,{xA:()=>c,yg:()=>d});var a=t(6540);function r(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function i(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function o(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?i(Object(t),!0).forEach((function(n){r(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):i(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,a,r=function(e,n){if(null==e)return{};var t,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)t=i[a],n.indexOf(t)>=0||(r[t]=e[t]);return r}(e,n);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)t=i[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(r[t]=e[t])}return r}var g=a.createContext({}),s=function(e){var n=a.useContext(g),t=n;return e&&(t="function"==typeof e?e(n):o(o({},n),e)),t},c=function(e){var n=s(e.components);return a.createElement(g.Provider,{value:n},e.children)},u="mdxType",p={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},m=a.forwardRef((function(e,n){var t=e.components,r=e.mdxType,i=e.originalType,g=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),u=s(t),m=r,d=u["".concat(g,".").concat(m)]||u[m]||p[m]||i;return t?a.createElement(d,o(o({ref:n},c),{},{components:t})):a.createElement(d,o({ref:n},c))}));function d(e,n){var t=arguments,r=n&&n.mdxType;if("string"==typeof e||r){var i=t.length,o=new Array(i);o[0]=m;var l={};for(var g in n)hasOwnProperty.call(n,g)&&(l[g]=n[g]);l.originalType=e,l[u]="string"==typeof e?e:r,o[1]=l;for(var s=2;s<i;s++)o[s]=t[s];return a.createElement.apply(null,o)}return a.createElement.apply(null,t)}m.displayName="MDXCreateElement"},4479:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>g,contentTitle:()=>o,default:()=>p,frontMatter:()=>i,metadata:()=>l,toc:()=>s});var a=t(8168),r=(t(6540),t(5680));t(1873);const i={},o="MSM",l={unversionedId:"icicle/golang-bindings/msm",id:"icicle/golang-bindings/msm",title:"MSM",description:"MSM Example",source:"@site/docs/icicle/golang-bindings/msm.md",sourceDirName:"icicle/golang-bindings",slug:"/icicle/golang-bindings/msm",permalink:"/icicle/golang-bindings/msm",editUrl:"https://github.com/ingonyama-zk/icicle/tree/main/docs/icicle/golang-bindings/msm.md",tags:[],version:"current",lastUpdatedBy:"VitaliiH",lastUpdatedAt:1715232285,formattedLastUpdatedAt:"5/9/2024",frontMatter:{},sidebar:"GettingStartedSidebar",previous:{title:"Golang bindings",permalink:"/icicle/golang-bindings"},next:{title:"MSM Pre computation",permalink:"/icicle/golang-bindings/msm-pre-computation"}},g={},s=[{value:"MSM Example",id:"msm-example",level:2},{value:"MSM Method",id:"msm-method",level:2},{value:"Parameters",id:"parameters",level:3},{value:"Return Value",id:"return-value",level:3},{value:"MSMConfig",id:"msmconfig",level:2},{value:"Fields",id:"fields",level:3},{value:"Default Configuration",id:"default-configuration",level:3},{value:"How do I toggle between the supported algorithms?",id:"how-do-i-toggle-between-the-supported-algorithms",level:2},{value:"How do I toggle between MSM modes?",id:"how-do-i-toggle-between-msm-modes",level:2},{value:"Support for G2 group",id:"support-for-g2-group",level:2}],c={toc:s},u="wrapper";function p(e){let{components:n,...t}=e;return(0,r.yg)(u,(0,a.A)({},c,t,{components:n,mdxType:"MDXLayout"}),(0,r.yg)("h1",{id:"msm"},"MSM"),(0,r.yg)("h2",{id:"msm-example"},"MSM Example"),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},'package main\n\nimport (\n  "github.com/ingonyama-zk/icicle/v2/wrappers/golang/core"\n  cr "github.com/ingonyama-zk/icicle/v2/wrappers/golang/cuda_runtime"\n  bn254 "github.com/ingonyama-zk/icicle/v2/wrappers/golang/curves/bn254"\n)\n\nfunc main() {\n  // Obtain the default MSM configuration.\n  cfg := bn254.GetDefaultMSMConfig()\n\n  // Define the size of the problem, here 2^18.\n  size := 1 << 18\n\n  // Generate scalars and points for the MSM operation.\n  scalars := bn254.GenerateScalars(size)\n  points := bn254.GenerateAffinePoints(size)\n\n  // Create a CUDA stream for asynchronous operations.\n  stream, _ := cr.CreateStream()\n  var p bn254.Projective\n\n  // Allocate memory on the device for the result of the MSM operation.\n  var out core.DeviceSlice\n  _, e := out.MallocAsync(p.Size(), p.Size(), stream)\n\n  if e != cr.CudaSuccess {\n    panic(e)\n  }\n\n  // Set the CUDA stream in the MSM configuration.\n  cfg.Ctx.Stream = &stream\n  cfg.IsAsync = true\n\n  // Perform the MSM operation.\n  e = bn254.Msm(scalars, points, &cfg, out)\n\n  if e != cr.CudaSuccess {\n    panic(e)\n  }\n\n  // Allocate host memory for the results and copy the results from the device.\n  outHost := make(core.HostSlice[bn254.Projective], 1)\n  cr.SynchronizeStream(&stream)\n  outHost.CopyFromDevice(&out)\n\n  // Free the device memory allocated for the results.\n  out.Free()\n}\n\n')),(0,r.yg)("h2",{id:"msm-method"},"MSM Method"),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},"func Msm(scalars core.HostOrDeviceSlice, points core.HostOrDeviceSlice, cfg *core.MSMConfig, results core.HostOrDeviceSlice) cr.CudaError\n")),(0,r.yg)("h3",{id:"parameters"},"Parameters"),(0,r.yg)("ul",null,(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"scalars")),": A slice containing the scalars for multiplication. It can reside either in host memory or device memory."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"points")),": A slice containing the points to be multiplied with scalars. Like scalars, these can also be in host or device memory."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"cfg")),": A pointer to an ",(0,r.yg)("inlineCode",{parentName:"li"},"MSMConfig")," object, which contains various configuration options for the MSM operation."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"results")),": A slice where the results of the MSM operation will be stored. This slice can be in host or device memory.")),(0,r.yg)("h3",{id:"return-value"},"Return Value"),(0,r.yg)("ul",null,(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"CudaError")),": Returns a CUDA error code indicating the success or failure of the MSM operation.")),(0,r.yg)("h2",{id:"msmconfig"},"MSMConfig"),(0,r.yg)("p",null,"The ",(0,r.yg)("inlineCode",{parentName:"p"},"MSMConfig")," structure holds configuration parameters for the MSM operation, allowing customization of its behavior to optimize performance based on the specifics of the operation or the underlying hardware."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},"type MSMConfig struct {\n    Ctx cr.DeviceContext\n    PrecomputeFactor int32\n    C int32\n    Bitsize int32\n    LargeBucketFactor int32\n    batchSize int32\n    areScalarsOnDevice bool\n    AreScalarsMontgomeryForm bool\n    arePointsOnDevice bool\n    ArePointsMontgomeryForm bool\n    areResultsOnDevice bool\n    IsBigTriangle bool\n    IsAsync bool\n}\n")),(0,r.yg)("h3",{id:"fields"},"Fields"),(0,r.yg)("ul",null,(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"Ctx")),": Device context containing details like device id and stream."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"PrecomputeFactor")),": Controls the number of extra points to pre-compute."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"C")),': Window bitsize, a key parameter in the "bucket method" for MSM.'),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"Bitsize")),": Number of bits of the largest scalar."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"LargeBucketFactor")),": Sensitivity to frequently occurring buckets."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"batchSize")),": Number of results to compute in one batch."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"areScalarsOnDevice")),": Indicates if scalars are located on the device."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"AreScalarsMontgomeryForm")),": True if scalars are in Montgomery form."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"arePointsOnDevice")),": Indicates if points are located on the device."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"ArePointsMontgomeryForm")),": True if point coordinates are in Montgomery form."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"areResultsOnDevice")),": Indicates if results are stored on the device."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"IsBigTriangle")),": If ",(0,r.yg)("inlineCode",{parentName:"li"},"true")," MSM will run in Large triangle accumulation if ",(0,r.yg)("inlineCode",{parentName:"li"},"false")," Bucket accumulation will be chosen. Default value: false."),(0,r.yg)("li",{parentName:"ul"},(0,r.yg)("strong",{parentName:"li"},(0,r.yg)("inlineCode",{parentName:"strong"},"IsAsync")),": If true, runs MSM asynchronously.")),(0,r.yg)("h3",{id:"default-configuration"},"Default Configuration"),(0,r.yg)("p",null,"Use ",(0,r.yg)("inlineCode",{parentName:"p"},"GetDefaultMSMConfig")," to obtain a default configuration, which can then be customized as needed."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},"func GetDefaultMSMConfig() MSMConfig\n")),(0,r.yg)("h2",{id:"how-do-i-toggle-between-the-supported-algorithms"},"How do I toggle between the supported algorithms?"),(0,r.yg)("p",null,"When creating your MSM Config you may state which algorithm you wish to use. ",(0,r.yg)("inlineCode",{parentName:"p"},"cfg.Ctx.IsBigTriangle = true")," will activate Large triangle accumulation and ",(0,r.yg)("inlineCode",{parentName:"p"},"cfg.Ctx.IsBigTriangle = false")," will activate Bucket accumulation."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},"...\n\n// Obtain the default MSM configuration.\ncfg := GetDefaultMSMConfig()\n\ncfg.Ctx.IsBigTriangle = true\n\n...\n")),(0,r.yg)("h2",{id:"how-do-i-toggle-between-msm-modes"},"How do I toggle between MSM modes?"),(0,r.yg)("p",null,"Toggling between MSM modes occurs automatically based on the number of results you are expecting from the ",(0,r.yg)("inlineCode",{parentName:"p"},"MSM")," function."),(0,r.yg)("p",null,"The number of results is interpreted from the size of ",(0,r.yg)("inlineCode",{parentName:"p"},"var out core.DeviceSlice"),". Thus its important when allocating memory for ",(0,r.yg)("inlineCode",{parentName:"p"},"var out core.DeviceSlice")," to make sure that you are allocating ",(0,r.yg)("inlineCode",{parentName:"p"},"<number of results> X <size of a single point>"),"."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},"... \n\nbatchSize := 3\nvar p G2Projective\nvar out core.DeviceSlice\nout.Malloc(batchSize*p.Size(), p.Size())\n\n...\n")),(0,r.yg)("h2",{id:"support-for-g2-group"},"Support for G2 group"),(0,r.yg)("p",null,"To activate G2 support first you must make sure you are building the static libraries with G2 feature enabled as described in the ",(0,r.yg)("a",{parentName:"p",href:"/icicle/golang-bindings#using-icicle-golang-bindings-in-your-project"},"Golang building instructions"),"."),(0,r.yg)("p",null,"Now you may import ",(0,r.yg)("inlineCode",{parentName:"p"},"g2")," package of the specified curve."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},'import (\n    "github.com/ingonyama-zk/icicle/v2/wrappers/golang/curves/bn254/g2"\n)\n')),(0,r.yg)("p",null,"This package include ",(0,r.yg)("inlineCode",{parentName:"p"},"G2Projective")," and ",(0,r.yg)("inlineCode",{parentName:"p"},"G2Affine")," points as well as a ",(0,r.yg)("inlineCode",{parentName:"p"},"G2Msm")," method."),(0,r.yg)("pre",null,(0,r.yg)("code",{parentName:"pre",className:"language-go"},'package main\n\nimport (\n  "github.com/ingonyama-zk/icicle/v2/wrappers/golang/core"\n  bn254 "github.com/ingonyama-zk/icicle/v2/wrappers/golang/curves/bn254"\n  g2 "github.com/ingonyama-zk/icicle/v2/wrappers/golang/curves/bn254/g2"\n)\n\nfunc main() {\n  cfg := bn254.GetDefaultMSMConfig()\n  size := 1 << 12\n  batchSize := 3\n  totalSize := size * batchSize\n  scalars := bn254.GenerateScalars(totalSize)\n  points := g2.G2GenerateAffinePoints(totalSize)\n\n  var p g2.G2Projective\n  var out core.DeviceSlice\n  out.Malloc(batchSize*p.Size(), p.Size())\n  g2.G2Msm(scalars, points, &cfg, out)\n}\n\n')),(0,r.yg)("p",null,(0,r.yg)("inlineCode",{parentName:"p"},"G2Msm")," works the same way as normal MSM, the difference is that it uses G2 Points."))}p.isMDXComponent=!0},1873:(e,n,t)=>{t(6540)}}]);