template <typename T> T schlick_fresnel_90(const T &F90, Real cos_theta) {
    return Real(1) + (F90 - Real(1)) * 
            pow(Real(1) - cos_theta, Real(5));
}

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    auto h_out = dot(h, dir_out);
    auto n_out_abs = fabs(dot(frame.n, dir_out));
    auto n_in_abs = fabs(dot(frame.n, dir_in));

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real F_D90 = .5 + 2 * roughness * h_out * h_out;
    Spectrum f_base_diffuse = base_color / c_PI * schlick_fresnel_90(F_D90, n_in_abs) * schlick_fresnel_90(F_D90, n_out_abs) * n_out_abs;
    
    Real F_SS90 = roughness * h_out * h_out;
    Spectrum f_subsurface = 1.25 * base_color / c_PI * (schlick_fresnel_90(F_SS90, n_in_abs) * schlick_fresnel_90(F_SS90, n_out_abs) *
                            (Real(1) / (n_in_abs + n_out_abs) - .5) + .5) * n_out_abs;

    return (Real(1) - subsurface) * f_base_diffuse + subsurface * f_subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fabs(dot(frame.n, dir_out)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
