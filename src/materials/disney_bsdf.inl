#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    // bool reflect = dot(vertex.geometric_normal, dir_in) *
    //                dot(vertex.geometric_normal, dir_out) > 0;
    bool reflect = dot(vertex.geometric_normal, dir_in) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;

    Spectrum C_tint = (luminance(base_color) > 0) ? (base_color / luminance(base_color)) : make_const_spectrum(1);
    Spectrum K_s = (1 - specular_tint) + specular_tint * C_tint;
    Real R_0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    Spectrum C_0 = specular * R_0 * (1 - metallic) * K_s + metallic * base_color;

    Spectrum f_diffuse, f_sheen, f_metal, f_clearcoat, f_glass;
    if (reflect) {
        f_diffuse = operator()(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface});
        f_sheen = operator()(DisneySheen{bsdf.base_color, bsdf.sheen_tint});
        // f_metal = eval_op::operator()(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic});
        f_clearcoat = operator()(DisneyClearcoat{bsdf.clearcoat_gloss});
        if (dot(vertex.geometric_normal, dir_in) < 0 ||
                dot(vertex.geometric_normal, dir_out) < 0) {
            // No light below the surface
            f_metal =  make_zero_spectrum();
        } else {
            // Homework 1: implement this!
            Vector3 h = normalize(dir_in + dir_out);
            auto h_out = dot(h, dir_out);
            auto n_in_abs = fabs(dot(frame.n, dir_in));

            Spectrum F_m = schlick_fresnel(C_0, h_out);

            Real aspect = sqrt(Real(1) - .9 * anisotropic);
            const Real alpha_min = 0.0001;
            Real alpha_x = fmax(alpha_min, roughness * roughness / aspect);
            Real a_x_2 = alpha_x * alpha_x;
            Real alpha_y = fmax(alpha_min, roughness * roughness * aspect);
            Real a_y_2 = alpha_y * alpha_y;
            Vector3 h_l = to_local(frame, h);
            Real D_m_idx = h_l.x * h_l.x / a_x_2 + h_l.y * h_l.y / a_y_2 + h_l.z * h_l.z;
            Real D_m = Real(1) / (c_PI * alpha_x * alpha_y * D_m_idx * D_m_idx);

            Vector3 omega_l_in = to_local(frame, dir_in);
            Vector3 omega_l_out = to_local(frame, dir_out);
            Real lambda_in = (sqrt(1 + (omega_l_in.x * omega_l_in.x * a_x_2 + omega_l_in.y * omega_l_in.y * a_y_2) /
                                (omega_l_in.z * omega_l_in.z)) - 1) / Real(2);
            Real lambda_out = (sqrt(1 + (omega_l_out.x * omega_l_out.x * a_x_2 + omega_l_out.y * omega_l_out.y * a_y_2) /
                                (omega_l_out.z * omega_l_out.z)) - 1) / Real(2);
            Real G_m = Real(1) / ((1 + lambda_in) * (1 + lambda_out));
            
            f_metal = F_m * D_m * G_m / (4 * n_in_abs);
        }
    } else {
        f_diffuse = make_zero_spectrum();
        f_sheen = make_zero_spectrum();
        f_metal = make_zero_spectrum();
        f_clearcoat = make_zero_spectrum();
    }
    f_glass = operator()(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta});
    
    return (1 - specular_transmission) * (1 - metallic) * f_diffuse +
            (1 - metallic) * sheen * f_sheen +
            (1 - specular_transmission * (1- metallic)) * f_metal + 
            .25 * clearcoat * f_clearcoat +
            (1 - metallic) * specular_transmission * f_glass;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // bool reflect = dot(vertex.geometric_normal, dir_in) *
    //                dot(vertex.geometric_normal, dir_out) > 0;
    bool reflect = dot(vertex.geometric_normal, dir_in) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real weight_diffuse = 0;
    Real weight_metal = 0;
    Real weight_glass = 0;
    Real weight_clearcoat = 0;
    if (reflect) {
        weight_diffuse = (1 - metallic) * (1 - specular_transmission);
        weight_metal = (1 - specular_transmission * (1 - metallic));
        weight_clearcoat = .25 * clearcoat;
    }
    weight_glass = (1 - metallic) * specular_transmission;
    Real weight_total = weight_diffuse + weight_metal + weight_glass + weight_clearcoat;
    if (weight_total < 0.00001) {
        return 0;
    }
    Real p_diffuse = operator()(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface});
    Real p_glass = operator()(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta});
    Real p_metal = operator()(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic});
    Real p_clearcoat = operator()(DisneyClearcoat{bsdf.clearcoat_gloss});
        
    
    return (p_diffuse * weight_diffuse + 
            p_glass * weight_glass +
            p_clearcoat * weight_clearcoat +
            p_metal * weight_metal) / weight_total;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    bool reflect = dot(vertex.geometric_normal, dir_in) > 0;

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real weight_diffuse = 0;
    Real weight_metal = 0;
    Real weight_glass = (1 - metallic) * specular_transmission;
    Real weight_clearcoat = 0;
    if (reflect) {
        weight_diffuse = (1 - metallic) * (1 - specular_transmission);
        weight_metal = (1 - specular_transmission * (1 - metallic));
        weight_clearcoat = .25 * clearcoat;
    }
    Real weight_total = weight_diffuse + weight_metal + weight_glass + weight_clearcoat;
    const Real threshold = 0.00001;
    if (weight_total < threshold) {
        return {};
    }
    Real rnd = rnd_param_w;
    Real select_glass = weight_glass / weight_total;
    Real select_diffuse = weight_diffuse / weight_total + select_glass;
    Real select_metal = weight_metal / weight_total + select_diffuse;
    // Real rnd = rnd_param_w * weight_total / weight_glass;
    // Real select_glass = Real(1);
    // Real select_diffuse = weight_diffuse / weight_glass + select_glass;
    // Real select_metal = weight_metal / weight_glass + select_diffuse;
    
    if (rnd <= select_glass) {
        // rnd_param_w = rnd * weight_total / weight_glass;
        // return operator()(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta});
        if (weight_glass >= threshold) {
            return sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                            dir_in, vertex, texture_pool, rnd_param_uv, rnd * weight_total / weight_glass, dir);
        }
    }
    if (rnd <= select_diffuse) {
        if (weight_diffuse >= threshold) {
            return operator()(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface});
        }
    } 
    if (rnd <= select_metal) {
        if (weight_metal >= threshold) {
            return operator()(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic});
        }
    }
    return operator()(DisneyClearcoat{bsdf.clearcoat_gloss});
    
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
