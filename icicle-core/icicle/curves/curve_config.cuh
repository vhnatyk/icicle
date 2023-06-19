#pragma once

#if defined(G2_DEFINED)
#include "../primitives/extension_field.cuh"
#else
#include "../primitives/field.cuh"
#endif
#include "../primitives/projective.cuh"

#include "bls12_381.cuh"
// #include "bn254.cuh"


typedef Field<fp_config> scalar_field_t;
typedef scalar_field_t scalar_t;
typedef Field<fq_config> point_field_t;
static constexpr point_field_t b = point_field_t{ weierstrass_b };
typedef Projective<point_field_t, scalar_field_t, b> projective_t;
typedef Affine<point_field_t> affine_t;
#if defined(G2_DEFINED)
typedef ExtensionField<fq_config> g2_point_field_t;
static constexpr g2_point_field_t b_g2 = g2_point_field_t{ point_field_t{ weierstrass_b_g2_re },  point_field_t{ weierstrass_b_g2_im }};
typedef Projective<g2_point_field_t, scalar_field_t, b_g2> g2_projective_t;
typedef Affine<g2_point_field_t> g2_affine_t;
#endif
