#include "../microfacet.h"
#include <cmath>
#define WAVELENGTH_R 630
#define WAVELENGTH_G 532
#define WAVELENGTH_B 465
Spectrum eval_op::operator()(const IridescentMicrofacet &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real filmEta = bsdf.filmEta;
    Real height = bsdf.height;

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real h_dot_in = dot(half_vector, dir_in);
    Vector3 v_par = normalize(dir_out - dir_in);
    auto n_in_abs = fmax(0, dot(frame.n, dir_in));

    Real sin = (1 - h_dot_in * h_dot_in) / filmEta / filmEta;
    Real cos = sqrt(1 - sin);
    sin = sqrt(sin);
    Vector3 dir_out_under = -cos * half_vector + sin * v_par;
    Real h_under_dot_in = -dot(half_vector, dir_out_under);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real R12 = fresnel_dielectric(h_dot_in, filmEta);
    Real T12 = 1 - R12;
    Real R21 = fresnel_dielectric(h_under_dot_in, 1.0 / filmEta);
    Real T21 = 1 - R21;
    Spectrum R23 = schlick_fresnel(base_color, h_under_dot_in);;
    Real delta = height * 2 * filmEta * cos;
    Real delta_phi_R = delta / WAVELENGTH_R * 2 * c_PI;
    Real delta_phi_G = delta / WAVELENGTH_G * 2 * c_PI;
    Real delta_phi_B = delta / WAVELENGTH_B * 2 * c_PI;
    Real cos_R = std::cos(delta_phi_R);
    Real cos_G = std::cos(delta_phi_G);
    Real cos_B = std::cos(delta_phi_B);
    Spectrum cos_I = Vector3(cos_R, cos_G, cos_B);
    Spectrum Y = T12 * R23 * T21;
    Spectrum Z = R21 * R23;
    Spectrum X = Y - R12 * Z;
    Spectrum R_I = (R12 * R12 + 2 * R12 * X * cos_I + X * X) / (1 - 2.0 * Z * cos_I + Z * Z);

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);
    Real D_m_idx = h_l.x * h_l.x / (alpha_x * alpha_x) + h_l.y * h_l.y / (alpha_y * alpha_y) + h_l.z * h_l.z;
    Real D_m = Real(1) / (c_PI * alpha_x * alpha_y * D_m_idx * D_m_idx);

    Vector3 omega_l_in = to_local(frame, dir_in);
    Vector3 omega_l_out = to_local(frame, dir_out);
    Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * alpha_x * alpha_x + omega_l_in.y * omega_l_in.y * alpha_y * alpha_y) /
                        (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
    Real lambda_out = (sqrt(1 + (omega_l_out.x * omega_l_out.x * alpha_x * alpha_x + omega_l_out.y * omega_l_out.y * alpha_y * alpha_y) /
                        (omega_l_out.z * omega_l_out.z)) - 1) / Real(2);
    Real G_m = Real(1) / ((1 + lambda_in) * (1 + lambda_out));
    
    return R_I * D_m * G_m / (4 * n_in_abs);
}

Real pdf_sample_bsdf_op::operator()(const IridescentMicrofacet &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);
    auto n_dot_in = dot(frame.n, dir_in);

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, h);
    Real D_m_idx = h_l.x * h_l.x / (alpha_x * alpha_x) + h_l.y * h_l.y / (alpha_y * alpha_y) + h_l.z * h_l.z;
    Real D_m = Real(1) / (c_PI * alpha_x * alpha_y * D_m_idx * D_m_idx);

    Vector3 omega_l_in = to_local(frame, dir_in);
    Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * alpha_x * alpha_x + omega_l_in.y * omega_l_in.y * alpha_y * alpha_y) /
                        (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
    Real G_m = Real(1) / (1 + lambda_in);
    Real spec_prob = (G_m * D_m) / (4 * n_dot_in);
    return spec_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const IridescentMicrofacet &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // We use the reflectance to choose between sampling the dielectric or diffuse layer.
    Spectrum Ks = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real lS = luminance(Ks);
    if (lS <= 0) {
        return {};
    }

    // Sample from the specular lobe.

    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Clamp roughness to avoid numerical issues.
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal =
        sample_visible_normals_anisotropic(local_dir_in, alpha_x, alpha_y, rnd_param_uv);
    
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const IridescentMicrofacet &bsdf) const {
    return bsdf.base_color;
}