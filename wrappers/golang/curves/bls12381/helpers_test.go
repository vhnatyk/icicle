package bls12381

import (
	"math/rand"
)

func generateRandomLimb(size int) []uint64 {
	limbs := make([]uint64, size)
	for i := range limbs {
		limbs[i] = rand.Uint64()
	}
	return limbs
}

func generateLimbOne(size int) []uint64 {
	limbs := make([]uint64, size)
	limbs[0] = 1
	return limbs
}

func generateBytesArray(size int) ([]byte, []uint64) {
	baseBytes := []byte{1, 2, 3, 4, 5, 6, 7, 8}
	var bytes []byte
	var limbs []uint64
	for i := 0; i < size; i++ {
		bytes = append(bytes, baseBytes...)
		limbs = append(limbs, 578437695752307201)
	}

	return bytes, limbs
}
