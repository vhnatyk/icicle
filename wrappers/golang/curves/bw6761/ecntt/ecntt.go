package ecntt

// #cgo CFLAGS: -I./include/
// #include "ecntt.h"
import "C"

import (
	"unsafe"

	"github.com/ingonyama-zk/icicle/wrappers/golang/core"
	cr "github.com/ingonyama-zk/icicle/wrappers/golang/cuda_runtime"
	icicle_bw6761 "github.com/ingonyama-zk/icicle/wrappers/golang/curves/bw6761"
)

func ECNtt[T any](points core.HostOrDeviceSlice, dir core.NTTDir, cfg *core.NTTConfig[T], results core.HostOrDeviceSlice) core.IcicleError {
	core.NttCheck[T](points, cfg, results)

	var pointsPointer unsafe.Pointer
	if points.IsOnDevice() {
		pointsPointer = points.(core.DeviceSlice).AsPointer()
	} else {
		pointsPointer = unsafe.Pointer(&points.(core.HostSlice[icicle_bw6761.Projective])[0])
	}
	cPoints := (*C.projective_t)(pointsPointer)
	cSize := (C.int)(points.Len() / int(cfg.BatchSize))
	cDir := (C.int)(dir)
	cCfg := (*C.NTTConfig)(unsafe.Pointer(cfg))

	var resultsPointer unsafe.Pointer
	if results.IsOnDevice() {
		resultsPointer = results.(core.DeviceSlice).AsPointer()
	} else {
		resultsPointer = unsafe.Pointer(&results.(core.HostSlice[icicle_bw6761.Projective])[0])
	}
	cResults := (*C.projective_t)(resultsPointer)

	__ret := C.bw6_761ECNTTCuda(cPoints, cSize, cDir, cCfg, cResults)
	err := (cr.CudaError)(__ret)
	return core.FromCudaError(err)
}
