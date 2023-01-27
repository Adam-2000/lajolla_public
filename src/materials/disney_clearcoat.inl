#include "../microfacet.h"
inline Real GTR1(Real n_dot_h, Real a) {
    if(a >= 1) {
        return c_INVPI;
    }
    Real a2 = a * a;
    return (a2 - Real(1)) / (c_PI * log2(a2) * (Real(1) + (a2 - Real(1)) * n_dot_h * n_dot_h));
}

// inline Real SeparableSmithGGXG1(Real n_dot_omega, float a) {
//     float a2 = a * a;
//     return Real(2) / (Real(1) + sqrt(a2 + (1 - a2) * n_dot_omega * n_dot_omega));
// }

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    auto n_in_abs = fabs(dot(frame.n, dir_in));

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    const Real eta = 1.5;
    Real r_0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    Real F_c = schlick_fresnel(r_0, h_out_abs);

    // Vector3 h_l = to_local(frame, h);
    Real alpha_g = fmax((Real(1) - clearcoat_gloss) * .1 + clearcoat_gloss * .001, 0.0001);
    Real D_c = GTR1(dot(frame.n, h), alpha_g);

    Vector3 omega_l_in = to_local(frame, dir_in);
    Vector3 omega_l_out = to_local(frame, dir_out);
    Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * .25 * .25 + omega_l_in.y * omega_l_in.y * .25 * .25) /
                        (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
    Real lambda_out = (sqrt(1 + (omega_l_out.x * omega_l_out.x * .25 * .25 + omega_l_out.y * omega_l_out.y * .25 * .25) /
                        (omega_l_out.z * omega_l_out.z)) - 1) / Real(2);
    Real G_c = Real(1) / ((1 + lambda_in) * (1 + lambda_out));
    // Real G_c = SeparableSmithGGXG1(omega_l_in.z, 0.25) * SeparableSmithGGXG1(omega_l_out.z, 0.25);
    
    return make_const_spectrum(F_c * D_c * G_c / (4 * n_in_abs));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real h_dot_out = dot(half_vector, dir_out);
    Real n_dot_out = dot(frame.n, dir_out);
    // Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    Real alpha = (Real(1) - clearcoat_gloss) * .1 + clearcoat_gloss * .001;
    alpha = std::clamp(alpha, Real(0.0001), Real(1));

    Real D = GTR1(n_dot_h, alpha);
    // (4 * cos_theta_v) is the Jacobian of the reflectiokn
    Real spec_prob = (D * n_dot_h) / (4 * h_dot_out);
    return spec_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (Real(1) - clearcoat_gloss) * .1 + clearcoat_gloss * .001;
    Real alpha_2 = alpha * alpha;

    Real cos_h_elevation = sqrt((1 - pow(alpha_2, 1 - rnd_param_uv.x)) / (1 - alpha_2));
    Real sin_h_elevation = sqrt(1 - cos_h_elevation * cos_h_elevation);
    Real h_azimuth = Real(2) * c_PI * rnd_param_uv.y;
    Real h_l_x = sin_h_elevation * cos(h_azimuth);
    Real h_l_y = sin_h_elevation * sin(h_azimuth);
    Real h_l_z = cos_h_elevation;

    Vector3 micro_normal = to_world(frame, Vector3(h_l_x, h_l_y, h_l_z));
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, micro_normal) * micro_normal);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, alpha /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
