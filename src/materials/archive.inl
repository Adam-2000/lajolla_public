
Spectrum eval_op::operator()(const ThinFilm &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real filmEta = bsdf.filmEta;
    Real height = bsdf.height;
    Vector3 half_vector = frame.n;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    Real h_dot_in = dot(half_vector, dir_in);
    auto n_in_abs = fmax(0, dot(frame.n, dir_in));

    Real sin = (1 - h_dot_in * h_dot_in) / filmEta / filmEta;
    Real cos = sqrt(1 - sin);
    sin = sqrt(sin);
    
    Real R = fresnel_dielectric(h_dot_in, filmEta);
    Real T = 1 - R;
    Real delta = height * 2 * filmEta * (2 / cos - cos);
    Real delta_phi_R = delta / WAVELENGTH_R * 2 * c_PI;
    Real delta_phi_G = delta / WAVELENGTH_G * 2 * c_PI;
    Real delta_phi_B = delta / WAVELENGTH_B * 2 * c_PI;
    Real cos_R = std::cos(delta_phi_R);
    Real cos_G = std::cos(delta_phi_G);
    Real cos_B = std::cos(delta_phi_B);
    Spectrum cos_I = Vector3(cos_R, cos_G, cos_B);
    Real Y = T * R * T;
    Real Z = R * R;
    Real X = Y - R * Z;
    Spectrum R_I = base_color * (R * R + 2 * R * X * cos_I + X * X) / (1 - 2.0 * Z * cos_I + Z * Z);
    if (reflect) {
        R_I = 1 - R_I;
    }
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

Real pdf_sample_bsdf_op::operator()(const ThinFilm &bsdf) const {
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

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = frame.n;
    }
    Real filmEta = bsdf.filmEta;
    Real height = bsdf.height;
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
   
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    Real h_dot_in = dot(half_vector, dir_in);

    Real sin = (1 - h_dot_in * h_dot_in) / filmEta / filmEta;
    Real cos = sqrt(1 - sin);
    sin = sqrt(sin);
    
    Real R = fresnel_dielectric(h_dot_in, filmEta);
    Real T = 1 - R;
    Real delta = height * 2 * filmEta * (2 / cos - cos);
    Real delta_phi_R = delta / WAVELENGTH_R * 2 * c_PI;
    Real delta_phi_G = delta / WAVELENGTH_G * 2 * c_PI;
    Real delta_phi_B = delta / WAVELENGTH_B * 2 * c_PI;
    Real cos_R = std::cos(delta_phi_R);
    Real cos_G = std::cos(delta_phi_G);
    Real cos_B = std::cos(delta_phi_B);
    Spectrum cos_I = Vector3(cos_R, cos_G, cos_B);
    Real Y = T * R * T;
    Real Z = R * R;
    Real X = Y - R * Z;
    Spectrum R_I = base_color * (R * R + 2 * R * X * cos_I + X * X) / (1 - 2.0 * Z * cos_I + Z * Z);
    if (reflect) {
        R_I = 1 - R_I;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

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

    return (average(R_I) * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const ThinFilm &bsdf) const {

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
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real filmEta = bsdf.filmEta;
    Real height = bsdf.height;
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
    Real h_dot_in = dot(half_vector, dir_in);

    Real sin = (1 - h_dot_in * h_dot_in) / filmEta / filmEta;
    Real cos = sqrt(1 - sin);
    sin = sqrt(sin);
    
    Real R = fresnel_dielectric(h_dot_in, filmEta);
    Real T = 1 - R;
    Real delta = height * 2 * filmEta * (2 / cos - cos);
    Real delta_phi_R = delta / WAVELENGTH_R * 2 * c_PI;
    Real delta_phi_G = delta / WAVELENGTH_G * 2 * c_PI;
    Real delta_phi_B = delta / WAVELENGTH_B * 2 * c_PI;
    Real cos_R = std::cos(delta_phi_R);
    Real cos_G = std::cos(delta_phi_G);
    Real cos_B = std::cos(delta_phi_B);
    Spectrum cos_I = Vector3(cos_R, cos_G, cos_B);
    Real Y = T * R * T;
    Real Z = R * R;
    Real X = Y - R * Z;
    Spectrum R_I = base_color * (R * R + 2 * R * X * cos_I + X * X) / (1 - 2.0 * Z * cos_I + Z * Z);
    if (rnd_param_w <= average(R_I)) {
        // Reflect over the world space normal
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, roughness /* roughness */
        };
    } else {
        return BSDFSampleRecord{
            -dir_in,
            Real(0) /* eta */, roughness /* roughness */
        };
    }
}

TextureSpectrum get_texture_op::operator()(const ThinFilm &bsdf) const {
    return bsdf.base_color;
}