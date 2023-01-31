#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Compute F / D / G
    // Note that we use the incoming direction
    // for evaluating the Fresnel reflection amount.
    // We can also use outgoing direction -- then we would need to
    // use 1/bsdf.eta and we will get the same result.
    // However, using the incoming direction allows
    // us to use F to decide whether to reflect or refract during sampling.
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    Real h_dot_in_abs = fabs(h_dot_in);
    Real h_dot_out_abs = fabs(h_dot_out);
    Real n_dot_in_abs = fabs(dot(frame.n, dir_in));
    // Real F_g = fresnel_dielectric(h_dot_in_abs, h_dot_out_abs, eta);
    Real F_g = fresnel_dielectric(h_dot_in, eta);
    // Real F_0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    // Real F_g = schlick_fresnel(F_0, dot(frame.n, dir_in));

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);
    Real D_g_idx = h_l.x * h_l.x / (alpha_x * alpha_x) + h_l.y * h_l.y / (alpha_y * alpha_y) + h_l.z * h_l.z;
    Real D_g = Real(1) / (c_PI * alpha_x * alpha_y * D_g_idx * D_g_idx);

    Vector3 omega_l_in = to_local(frame, dir_in);
    Vector3 omega_l_out = to_local(frame, dir_out);
    Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * alpha_x * alpha_x + omega_l_in.y * omega_l_in.y * alpha_y * alpha_y) /
                        (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
    Real lambda_out = (sqrt(1 + (omega_l_out.x * omega_l_out.x * alpha_x * alpha_x + omega_l_out.y * omega_l_out.y * alpha_y * alpha_y) /
                        (omega_l_out.z * omega_l_out.z)) - 1) / Real(2);
    Real G_g = Real(1) / ((1 + lambda_in) * (1 + lambda_out));
    if (reflect) {
        return base_color * (F_g * D_g * G_g) / (4 * n_dot_in_abs);
    } else {
        // Snell-Descartes law predicts that the light will contract/expand 
        // due to the different index of refraction. So the normal BSDF needs
        // to scale with 1/eta^2. However, the "adjoint" of the BSDF does not have
        // the eta term. This is due to the non-reciprocal nature of the index of refraction:
        // f(wi -> wo) / eta_o^2 = f(wo -> wi) / eta_i^2
        // thus f(wi -> wo) = f(wo -> wi) (eta_o / eta_i)^2
        // The adjoint of a BSDF is defined as swapping the parameter, and
        // this cancels out the eta term.
        // See Chapter 5 of Eric Veach's thesis "Robust Monte Carlo Methods for Light Transport Simulation"
        // for more details.
        return sqrt(base_color) * (1 - F_g) * D_g * G_g * h_dot_in_abs * h_dot_out_abs /
                (n_dot_in_abs * (h_dot_in + eta * h_dot_out) * (h_dot_in + eta * h_dot_out));
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    // Real F_0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    // Real F = schlick_fresnel(F_0, dot(frame.n, dir_in));
    // auto n_dot_in = dot(frame.n, dir_in);

    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 h_l = to_local(frame, half_vector);
    Real D_m_idx = h_l.x * h_l.x / (alpha_x * alpha_x) + h_l.y * h_l.y / (alpha_y * alpha_y) + h_l.z * h_l.z;
    Real D = Real(1) / (c_PI * alpha_x * alpha_y * D_m_idx * D_m_idx);

    Vector3 omega_l_in = to_local(frame, dir_in);
    Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * alpha_x * alpha_x + omega_l_in.y * omega_l_in.y * alpha_y * alpha_y) /
                        (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
    Real G_in = Real(1) / (1 + lambda_in);

    if (reflect) {
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(Real(1) - .9 * anisotropic);
    const Real alpha_min = 0.0001;
    Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
    Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal =
        sample_visible_normals_anisotropic(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    // Real F_0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    // Real F = schlick_fresnel(F_0, dot(frame.n, dir_in));

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
