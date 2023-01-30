#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);
    auto h_out_abs = fabs(dot(h, dir_out));
    auto n_in_abs = fmax(0, dot(frame.n, dir_in));

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Spectrum F_m = schlick_fresnel(base_color, h_out_abs);

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, h);
    Real D_m_idx = h_l.x * h_l.x / (alpha_x * alpha_x) + h_l.y * h_l.y / (alpha_y * alpha_y) + h_l.z * h_l.z;
    Real D_m = Real(1) / (c_PI * alpha_x * alpha_y * D_m_idx * D_m_idx);

    Vector3 omega_l_in = to_local(frame, dir_in);
    Vector3 omega_l_out = to_local(frame, dir_out);
    Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * alpha_x * alpha_x + omega_l_in.y * omega_l_in.y * alpha_y * alpha_y) /
                        (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
    Real lambda_out = (sqrt(1 + (omega_l_out.x * omega_l_out.x * alpha_x * alpha_x + omega_l_out.y * omega_l_out.y * alpha_y * alpha_y) /
                        (omega_l_out.z * omega_l_out.z)) - 1) / Real(2);
    Real G_m = Real(1) / ((1 + lambda_in) * (1 + lambda_out));
    
    return F_m * D_m * G_m / (4 * n_in_abs);
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    Spectrum S = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real lS = luminance(S);
    if (lS <= 0) {
        return 0;
    }
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // We use the reflectance to determine whether to choose specular sampling lobe or diffuse.

    // For the specular lobe, we use the ellipsoidal sampling from Heitz 2018
    // "Sampling the GGX Distribution of Visible Normals"
    // https://jcgt.org/published/0007/04/01/
    // this importance samples smith_masking(cos_theta_in) * GTR2(cos_theta_h, roughness) * cos_theta_out
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness);
    Real D = GTR2(n_dot_h, roughness);
    // (4 * cos_theta_v) is the Jacobian of the reflectiokn
    Real spec_prob = (G * D) / (4 * n_dot_in);
    return spec_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
