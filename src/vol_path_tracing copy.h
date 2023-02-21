#pragma once
#include "spectrum.h"
#include "scene.h"
#include "pcg.h"
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng);
// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        // Hit background. Account for the environment map if needed.
        if (has_envmap(scene)) {
            const Light &envmap = get_envmap(scene);
            return emission(envmap,
                            -ray.dir, // pointing outwards from light
                            ray_diff.spread,
                            PointAndNormal{}, // dummy parameter for envmap
                            scene);
        }
        return make_zero_spectrum();
    }
    PathVertex vertex = *vertex_;

    Spectrum radiance = make_zero_spectrum();
    int medium_id = dot(ray.dir, vertex.geometric_normal) > 0 ? vertex.interior_medium_id : vertex.exterior_medium_id;
    Spectrum sigma_a = get_sigma_a(scene.media[medium_id], vertex.position);
    Real t = dot(vertex.position - ray.org, ray.dir);
    Spectrum transmittance = exp(-sigma_a * t);

    // We hit a light immediately. 
    // This path has only two vertices and has contribution
    // C = W(v0, v1) * G(v0, v1) * L(v0, v1)
    if (is_light(scene.shapes[vertex.shape_id])) {
        radiance += emission(vertex, -ray.dir, scene);
    }
    return radiance * transmittance;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    Real u = next_pcg32_real<Real>(rng);
    int medium_id;
    if (!vertex_) {
        medium_id = scene.camera.medium_id;
    } else {
        medium_id = dot(ray.dir, vertex_->geometric_normal) > 0 ? vertex_->interior_medium_id : vertex_->exterior_medium_id;
    }
    Medium medium = scene.media[medium_id]; 
    Real sigma_t = get_majorant(medium, ray).x;
    Real t = -log(1 - u) / sigma_t;
    Real t_hit = vertex_.has_value() ? dot(vertex_->position - ray.org, ray.dir) : MAXFLOAT;
    Spectrum radiance = make_zero_spectrum();
    Vector3 dir_view = -ray.dir;
    Real transmittance, trans_pdf;
    if (t >= t_hit) {
    // if (vertex_.has_value()) {
        PathVertex vertex = *vertex_;
        // transmittance = exp(-sigma_t * t_hit);
        // trans_pdf = exp(-sigma_t * t_hit);
        if (is_light(scene.shapes[vertex.shape_id])) {
            radiance += emission(vertex, dir_view, scene);
        }
        // return radiance * transmittance / trans_pdf;
        return radiance;
    }

    transmittance = exp(-sigma_t * t);
    trans_pdf = transmittance * sigma_t;
    Vector3 p = ray.org + t * ray.dir;

    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
    PhaseFunction rho = get_phase_function(medium);
    
    Real G = 0;
    Real dist = distance(point_on_light.position, p);
    Vector3 dir_light;
    // if (!is_envmap(light)) {
        dir_light = normalize(point_on_light.position - p);
        Ray shadow_ray{p, dir_light, 
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                            dist};
        if (!occluded(scene, shadow_ray)) {
            G = abs(dot(dir_light, point_on_light.normal)) /
                distance_squared(point_on_light.position, p);
        }

    Real L_s1_pdf = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, p, scene);
    if (G > 0 && L_s1_pdf > 0) {
        Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Spectrum rho_value = eval(rho, dir_view, dir_light);
        Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        Spectrum L_s1_estimate = rho_value * L * exp(-sigma_t * dist) * G;
        radiance = (transmittance / trans_pdf) * get_sigma_s(medium, p) * L_s1_estimate / L_s1_pdf;
    }
    return radiance;
}

void update_medium_id(const PathVertex& vertex, const Ray& ray, int& medium_id) {
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        if (dot(ray.dir, vertex.geometric_normal) > 0) {
            medium_id = vertex.exterior_medium_id;
        } else {
            medium_id = vertex.interior_medium_id;
        }
    }
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    return vol_path_tracing_5(scene, x, y, rng);
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    // PathVertex vertex = *vertex_;
    // int curr_medium_id = scene.camera.medium_id;
    int curr_medium_id = scene.camera.medium_id; 
    Spectrum radiance = make_zero_spectrum();
    int bounces = -1;
    Spectrum current_path_throughput = make_const_spectrum(1);
    int max_depth = scene.options.max_depth;
    Vector3 next_org = ray.org;
    Vector3 next_dir = ray.dir;
    while (true) {
        ray.org = next_org;
        ray.dir = next_dir;
        bounces += 1;
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = vertex_.has_value() ? dot(vertex_->position - ray.org, ray.dir) : MAXFLOAT;
        Real transmittance = 1;
        Real trans_pdf = 1;
        if (curr_medium_id != -1) {
            Medium medium = scene.media[curr_medium_id];
            Real sigma_t = get_majorant(medium, ray).x;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = transmittance * sigma_t;
                next_org = ray.org + t * ray.dir;
            } else {
                transmittance = exp(-sigma_t * t_hit);
                trans_pdf = transmittance;
                next_org = ray.org + t_hit * ray.dir;
            }
        }
        current_path_throughput *= transmittance / trans_pdf;
        Vector3 dir_view = -ray.dir;
        if (!scatter && vertex_ && is_light(scene.shapes[vertex_->shape_id])) {
            radiance += current_path_throughput * emission(*vertex_, dir_view, scene);
        }
        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }
        if (!scatter && vertex_) {
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            }
        }
        if (scatter) {
            Medium medium = scene.media[curr_medium_id];
            PhaseFunction rho = get_phase_function(medium);
            Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(rho, dir_view, rho_rnd_param_uv);
            if (next_dir_) {
                next_dir = *next_dir_;
                current_path_throughput *= 
                    eval(rho, dir_view, next_dir) /
                    pdf_sample_phase(rho, dir_view, next_dir) *
                    get_sigma_s(medium, ray.org);
            } 
        } else {
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
    }
    return radiance;
}
// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    return vol_path_tracing_5(scene, x, y, rng);
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    // PathVertex vertex = *vertex_;
    // int curr_medium_id = scene.camera.medium_id;
    int curr_medium_id = scene.camera.medium_id; 
    Spectrum radiance = make_zero_spectrum();
    int bounces = -1;
    Spectrum current_path_throughput = make_const_spectrum(1);
    int max_depth = scene.options.max_depth;
    Vector3 next_org = ray.org;
    Vector3 next_dir = ray.dir;
    bool never_scatter = true;
    Real dir_pdf = 1;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1;
    while (true) {
        ray.org = next_org;
        ray.dir = next_dir;
        bounces += 1;
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = vertex_.has_value() ? dot(vertex_->position - ray.org, ray.dir) : MAXFLOAT;
        Real transmittance = 1;
        Real trans_pdf = 1;
        Real pdf_nee;
        Real G;
        Vector3 dir_view = -ray.dir;
        if (curr_medium_id == -1) {
            if (!vertex_.has_value()) {
                break;
            }
            next_org = ray.org + (t_hit + get_intersection_epsilon(scene)) * ray.dir;
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            }
            radiance += current_path_throughput * emission(*vertex_, dir_view, scene);
            break;
        }
        Medium medium = scene.media[curr_medium_id];
        Real sigma_t = get_majorant(medium, ray).x;
        Real u = next_pcg32_real<Real>(rng);
        Real t = -log(1 - u) / sigma_t;
        if (t < t_hit) {
            scatter = true;
            transmittance = exp(-sigma_t * t);
            trans_pdf = transmittance * sigma_t;
            next_org = ray.org + t * ray.dir;
        } else {
            transmittance = exp(-sigma_t * t_hit);
            trans_pdf = transmittance;
            next_org = ray.org + (t_hit + get_intersection_epsilon(scene)) * ray.dir;
        }
        current_path_throughput *= transmittance / trans_pdf;
        multi_trans_pdf *= trans_pdf;
        nee_p_cache = next_org;
        /* Sample Light */
        {
        Vector3 p = nee_p_cache;
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        // int light_id = 1;
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
        Vector3 p_prime = point_on_light.position;
        Real T_light = 1;
        int shadow_medium_id = curr_medium_id;
        int shadow_bounces = 0;
        Real p_trans_dir = 1;
        Vector3 dir_light = normalize(p_prime - p);
        Vector3 curr_p = p;
        while (true) {
            Real dist = distance(p_prime, curr_p);
            Ray shadow_ray{curr_p, dir_light, 
                            get_shadow_epsilon(scene),
                            (1 - get_shadow_epsilon(scene)) *
                                dist};
            std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray);
            Real next_t = dist;
            if (shadow_vertex_) {
                next_t = distance(curr_p, shadow_vertex_->position);
            }
            if (shadow_medium_id != -1) {
                Real shadow_sigma_t = get_majorant(scene.media[shadow_medium_id], shadow_ray).x;
                T_light *= exp(-shadow_sigma_t * next_t);
                p_trans_dir *= exp(-shadow_sigma_t * next_t);
            }
            if (!shadow_vertex_) {
                break;
            }
            if (shadow_vertex_->material_id >= 0) {
                T_light = -1;
                break;
            }
            shadow_bounces += 1;
            if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                T_light = -1;
                break;
            }
            update_medium_id(*shadow_vertex_, shadow_ray, shadow_medium_id);
            curr_p += next_t * dir_light;
        }
        if (T_light > 0) {
            PhaseFunction rho = get_phase_function(medium);
            Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Spectrum rho_value = eval(rho, dir_view, dir_light);
            Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
            pdf_nee = light_pmf(scene, light_id) *
                        pdf_point_on_light(light, point_on_light, p, scene);
            G = abs(dot(dir_light, point_on_light.normal)) /
                    distance_squared(p_prime, p);
            Spectrum contrib = T_light * G * rho_value * L / pdf_nee;
            // Spectrum contrib = T_light * rho_value * L / pdf_nee;
            // Spectrum contrib = T_light * rho_value * L / pdf_nee;
            Real pdf_phase = pdf_sample_phase(rho, dir_view, dir_light) * G * p_trans_dir;
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            // radiance += current_path_throughput * get_sigma_s(medium, p) * w * contrib;
            radiance += current_path_throughput * get_sigma_s(medium, p) * w * contrib;
        }
        }

        /* End Sampling Light */

        if (!scatter && vertex_ && is_light(scene.shapes[vertex_->shape_id])) {
            if (never_scatter) {
                radiance += current_path_throughput * emission(*vertex_, dir_view, scene);
            } else {
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                radiance += current_path_throughput * w * emission(*vertex_, dir_view, scene);
            }
            break;
        }
        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }
        if (!scatter && vertex_) {
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            } 
        }
        if (scatter) {
            never_scatter = false;
            PhaseFunction rho = get_phase_function(medium);
            Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(rho, dir_view, rho_rnd_param_uv);
            if (next_dir_) {
                next_dir = *next_dir_;
                dir_pdf = pdf_sample_phase(rho, dir_view, next_dir);
                multi_trans_pdf = 1;
                current_path_throughput *= 
                    eval(rho, dir_view, next_dir) /
                    dir_pdf *
                    get_sigma_s(medium, ray.org);
            }
        } else {
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
    }
    return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    // RayDifferential ray_diff = init_ray_differential(w, h);

    // PathVertex vertex = *vertex_;
    // int curr_medium_id = scene.camera.medium_id;
    int curr_medium_id = scene.camera.medium_id; 
    Spectrum radiance = make_zero_spectrum();
    int bounces = -1;
    Spectrum current_path_throughput = make_const_spectrum(1);
    int max_depth = scene.options.max_depth;
    Vector3 next_org = ray.org;
    Vector3 next_dir = ray.dir;
    bool never_scatter = true;
    Real dir_pdf = 1;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1;
    // const int xx = 312;
    // const int yy = 256;
    ray.tnear = get_intersection_epsilon(scene);
    while (true) {
        ray.org = next_org;
        ray.dir = next_dir;
        bounces += 1;
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);
        Real t_hit = vertex_.has_value() ? dot(vertex_->position - ray.org, ray.dir) : MAXFLOAT;
        Real transmittance = 1;
        Real trans_pdf = 1;
        // Real pdf_nee;
        Vector3 dir_view = -ray.dir;
        if (curr_medium_id == -1) {
            if (!vertex_.has_value()) {
                break;
            }
            next_org = ray.org + t_hit * ray.dir;
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            }
        } else {
            Medium medium = scene.media[curr_medium_id];
            Real sigma_t = get_majorant(medium, ray).x;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1 - u) / sigma_t;
            if (t < t_hit) {
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = transmittance * sigma_t;
                next_org = ray.org + t * ray.dir;
            } else {
                transmittance = exp(-sigma_t * t_hit);
                trans_pdf = transmittance;
                next_org = ray.org + t_hit * ray.dir;
            }
            current_path_throughput *= transmittance / trans_pdf;
            multi_trans_pdf *= trans_pdf;
        }
        /* Sample Light */
        if (scatter || vertex_->material_id != -1) {
            nee_p_cache = next_org;
            Vector3 p = nee_p_cache;
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            // int light_id = 1;
            const Light &light = scene.lights[light_id];
            PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
            Vector3 p_prime = point_on_light.position;
            Real T_light = 1;
            int shadow_medium_id = curr_medium_id;
            int shadow_bounces = 0;
            Real p_trans_dir = 1;
            Vector3 dir_light = normalize(p_prime - p);
            Vector3 curr_p = p;
            while (true) {
                Real dist = distance(p_prime, curr_p);
                Ray shadow_ray{curr_p, dir_light, 
                                10 * get_shadow_epsilon(scene),
                                // get_shadow_epsilon(scene),
                                (1 - get_shadow_epsilon(scene)) *
                                    dist};
                std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray);
                Real next_t = dist;
                if (shadow_vertex_) {
                    next_t = distance(curr_p, shadow_vertex_->position);
                }
                if (shadow_medium_id != -1) {
                    Real shadow_sigma_t = get_majorant(scene.media[shadow_medium_id], shadow_ray).x;
                    T_light *= exp(-shadow_sigma_t * next_t);
                    p_trans_dir *= exp(-shadow_sigma_t * next_t);
                }
                if (!shadow_vertex_) {
                    break;
                }
                if (shadow_vertex_->material_id >= 0) {
                    T_light = -1;
                    break;
                }
                shadow_bounces += 1;
                if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                    T_light = -1;
                    break;
                }
                update_medium_id(*shadow_vertex_, shadow_ray, shadow_medium_id);
                curr_p += next_t * dir_light;
            }
            if (T_light > 0) {
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
                Real pdf_nee = light_pmf(scene, light_id) *
                            pdf_point_on_light(light, point_on_light, p, scene);
                Real G_shadow = abs(dot(dir_light, point_on_light.normal)) /
                        distance_squared(p_prime, p);
                Spectrum f;
                Real pdf_2;
                Spectrum ratio = make_const_spectrum(1);
                if (scatter) {
                    Medium medium = scene.media[curr_medium_id];
                    PhaseFunction rho = get_phase_function(medium);
                    Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                    f = eval(rho, dir_view, dir_light);
                    pdf_2 = pdf_sample_phase(rho, dir_view, dir_light) * G_shadow * p_trans_dir;
                    ratio = get_sigma_s(medium, p);
                } else {
                    Material mat = scene.materials[vertex_->material_id];
                    f = eval(mat, dir_view, dir_light, *vertex_, scene.texture_pool);
                    pdf_2 = pdf_sample_bsdf(mat, dir_view, dir_light, *vertex_, scene.texture_pool) * G_shadow;
                }
                
                Spectrum contrib = T_light * G_shadow * f * L / pdf_nee;
                Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_2 * pdf_2);
                // Real w = 1;
                radiance += current_path_throughput * ratio * w * contrib;
            }
        }

        /* End Sampling Light */

        if (!scatter && is_light(scene.shapes[vertex_->shape_id])) {
            Spectrum L = emission(*vertex_, dir_view, scene);
            Real w = 1;
            if (!never_scatter) {
                Real G = abs(dot(ray.dir, vertex_->geometric_normal)) /
                        distance_squared(ray.org, vertex_->position);
                int light_id = get_area_light_id(scene.shapes[vertex_->shape_id]);
                assert(light_id >= 0);
                const Light &light = scene.lights[light_id];
                PointAndNormal light_point{vertex_->position, vertex_->geometric_normal};
                Real pdf_nee = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, light_point, nee_p_cache, scene);
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                // w = 1;
            }
            radiance += current_path_throughput * w * L;
            break;
        }
        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }
        if (!scatter) {
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            } 
            never_scatter = false;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            Material mat = scene.materials[vertex_->material_id];
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            *vertex_,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            next_dir = bsdf_sample_->dir_out;
            if (dot(vertex_->geometric_normal, next_dir) * dot(vertex_->geometric_normal, ray.dir) > 0) {
                update_medium_id(*vertex_, ray, curr_medium_id);
            }
            Spectrum f = eval(mat, dir_view, next_dir, *vertex_, scene.texture_pool);
            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, next_dir, *vertex_, scene.texture_pool);
            current_path_throughput *= f / bsdf_pdf;
            multi_trans_pdf = 1;
            dir_pdf = bsdf_pdf;
        } else {
            never_scatter = false;
            Medium medium = scene.media[curr_medium_id];
            PhaseFunction rho = get_phase_function(medium);
            Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(rho, dir_view, rho_rnd_param_uv);
            if (next_dir_.has_value()) {
                next_dir = next_dir_.value();
                dir_pdf = pdf_sample_phase(rho, dir_view, next_dir);
                multi_trans_pdf = 1;
                current_path_throughput *= 
                    eval(rho, dir_view, next_dir) /
                    dir_pdf *
                    get_sigma_s(medium, ray.org);
            }
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
    }
    return radiance;
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    // RayDifferential ray_diff = init_ray_differential(w, h);

    // PathVertex vertex = *vertex_;
    // int curr_medium_id = scene.camera.medium_id;
    int curr_medium_id = scene.camera.medium_id; 
    Spectrum radiance = make_zero_spectrum();
    int bounces = -1;
    Spectrum current_path_throughput = make_const_spectrum(1);
    int max_depth = scene.options.max_depth;
    Vector3 next_org = ray.org;
    Vector3 next_dir = ray.dir;
    bool never_scatter = true;
    Real dir_pdf = 1;
    Vector3 nee_p_cache;
    Vector3 nee_pdf_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    const int xx = 120;
    const int yy = 120;
    ray.tnear = get_intersection_epsilon(scene);
    int max_null_collisions = scene.options.max_null_collisions;
    if (x == xx && y == yy) {
        std::cout<<"start: " << std::endl;
    }
    while (true) {
        if (x == xx && y == yy) {
            std::cout<< "radiance: " << radiance << std::endl;
        }
        ray.org = next_org;
        ray.dir = next_dir;
        bounces += 1;
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);
        Real t_hit = vertex_.has_value() ? dot(vertex_->position - ray.org, ray.dir) : MAXFLOAT;
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        Spectrum trans_nee_pdf = make_const_spectrum(1);
        // Real pdf_nee;
        Vector3 dir_view = -ray.dir;
        if (curr_medium_id == -1) {
            if (!vertex_.has_value()) {
                break;
            }
            next_org = ray.org + t_hit * ray.dir;
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            }
        } else {
            Medium medium = scene.media[curr_medium_id];
            Spectrum majorant = get_majorant(medium, ray);
            Real u = next_pcg32_real<Real>(rng);
            int channel = max(min(int(floor(u * 3)), 2), 0);
            Real accum_t = 0;
            int iteration = 0;
            while (majorant[channel] > 0 && iteration++ < max_null_collisions) {
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;
                accum_t = min(accum_t + t, t_hit);
                Spectrum sigma_t = get_sigma_a(medium, ray.org + ray.dir * accum_t) + get_sigma_s(medium, ray.org + ray.dir * accum_t);
                Spectrum sigma_n = majorant - sigma_t;
                Spectrum exp_minus_m_t = exp(-majorant * t);
                Real max_m = max(majorant);
                if (t < dt) {
                    Spectrum real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        scatter = true;
                        transmittance *= exp_minus_m_t / max_m;
                        trans_dir_pdf *= exp_minus_m_t * majorant * real_prob / max_m;
                        break;
                    } else {
                        transmittance *= exp_minus_m_t * sigma_n / max_m;
                        trans_dir_pdf *= exp_minus_m_t * majorant * (1 - real_prob) / max_m;
                        trans_nee_pdf *= exp_minus_m_t * majorant / max_m;
                    }
                } else {
                    transmittance *= exp_minus_m_t;
                    trans_dir_pdf *= exp_minus_m_t;
                    trans_nee_pdf *= exp_minus_m_t;
                    break;
                }
            }
            next_org = ray.org + accum_t * ray.dir;
            // current_path_throughput *= transmittance / average(trans_dir_pdf);
            current_path_throughput *= transmittance / trans_dir_pdf;
            multi_trans_pdf *= trans_dir_pdf;
        }
        /* Sample Light */
        if (scatter || vertex_->material_id != -1) {
            nee_p_cache = next_org;
            nee_pdf_cache = trans_nee_pdf;
            Vector3 p = nee_p_cache;
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            // int light_id = 1;
            const Light &light = scene.lights[light_id];
            PointAndNormal point_on_light = sample_point_on_light(light, p, light_uv, shape_w, scene);
            Vector3 p_prime = point_on_light.position;
            Spectrum T_light = make_const_spectrum(1);
            int shadow_medium_id = curr_medium_id;
            int shadow_bounces = 0;
            Spectrum p_trans_nee = make_const_spectrum(1);
            Spectrum p_trans_dir = make_const_spectrum(1);
            Vector3 dir_light = normalize(p_prime - p);
            Vector3 curr_p = p;
            while (true) {
                Real dist = distance(p_prime, curr_p);
                Ray shadow_ray{curr_p, dir_light, 
                                10 * get_shadow_epsilon(scene),
                                // get_shadow_epsilon(scene),
                                (1 - get_shadow_epsilon(scene)) *
                                    dist};
                std::optional<PathVertex> shadow_vertex_ = intersect(scene, shadow_ray);
                Real next_t = dist;
                if (shadow_vertex_) {
                    next_t = distance(curr_p, shadow_vertex_->position);
                }
                if (shadow_medium_id != -1) {
                    Spectrum majorant = get_majorant(scene.media[shadow_medium_id], shadow_ray);
                    Medium medium = scene.media[shadow_medium_id];
                    Real u = next_pcg32_real<Real>(rng);
                    Real channel = max(min(int(floor(u * 3)), 2), 0);
                    Real accum_t = 0;
                    int iteration = 0;
                    while (majorant[channel] > 0 && iteration++ < max_null_collisions) {
                        Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                        Real dt = next_t - accum_t;
                        accum_t = min(accum_t + t, next_t);
                        Spectrum sigma_t = get_sigma_a(medium, shadow_ray.org + shadow_ray.dir * accum_t) + get_sigma_s(medium, shadow_ray.org + shadow_ray.dir * accum_t);
                        Spectrum sigma_n = majorant - sigma_t;
                        Spectrum exp_minus_m_t = exp(-majorant * t);
                        Real max_m = max(majorant);
                        if (t < dt) {
                            T_light *= exp_minus_m_t * sigma_n / max_m;
                            p_trans_nee *= exp_minus_m_t * majorant / max_m;
                            Spectrum real_prob = sigma_t / majorant;
                            p_trans_dir *= exp_minus_m_t * majorant * (1 - real_prob) / max_m;
                            if (max(T_light) < 0) {
                                break;
                            }
                        } else {
                            T_light *= exp_minus_m_t;
                            p_trans_dir *= exp_minus_m_t;
                            p_trans_nee *= exp_minus_m_t;
                            break;
                        }
                    }
                }
                if (!shadow_vertex_) {
                    break;
                }
                if (shadow_vertex_->material_id >= 0) {
                    T_light = make_const_spectrum(-1);
                    break;
                }
                shadow_bounces += 1;
                if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                    T_light = make_const_spectrum(-1);
                    break;
                }
                update_medium_id(*shadow_vertex_, shadow_ray, shadow_medium_id);
                curr_p += next_t * dir_light;
            }
            if (min(T_light) >= 0 && min(p_trans_nee) >= 0) {
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
                Spectrum pdf_nee = light_pmf(scene, light_id) *
                            pdf_point_on_light(light, point_on_light, p, scene)
                            * p_trans_nee;
                Real G_shadow = abs(dot(dir_light, point_on_light.normal)) /
                        distance_squared(p_prime, p);
                Spectrum f;
                Spectrum pdf_2;
                Spectrum ratio = make_const_spectrum(1);
                if (scatter) {
                    Medium medium = scene.media[curr_medium_id];
                    PhaseFunction rho = get_phase_function(medium);
                    Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                    f = eval(rho, dir_view, dir_light);
                    pdf_2 = pdf_sample_phase(rho, dir_view, dir_light) * G_shadow * p_trans_dir;
                    ratio = get_sigma_s(medium, p);
                } else {
                    Material mat = scene.materials[vertex_->material_id];
                    f = eval(mat, dir_view, dir_light, *vertex_, scene.texture_pool);
                    pdf_2 = make_const_spectrum(pdf_sample_bsdf(mat, dir_view, dir_light, *vertex_, scene.texture_pool) * G_shadow) * p_trans_dir;
                }
                
                // Spectrum contrib = T_light * G_shadow * f * L / average(pdf_nee);
                Spectrum contrib = T_light * G_shadow * f * L / pdf_nee;
                Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_2 * pdf_2);
                // Spectrum w = make_const_spectrum(1);
                radiance += current_path_throughput * ratio * w * contrib;
            }
        }

        /* End Sampling Light */

        if (!scatter && is_light(scene.shapes[vertex_->shape_id])) {
            Spectrum L = emission(*vertex_, dir_view, scene);
            Spectrum w = make_const_spectrum(1);
            if (!never_scatter) {
                Real G = abs(dot(ray.dir, vertex_->geometric_normal)) /
                        distance_squared(ray.org, vertex_->position);
                int light_id = get_area_light_id(scene.shapes[vertex_->shape_id]);
                assert(light_id >= 0);
                const Light &light = scene.lights[light_id];
                PointAndNormal light_point{vertex_->position, vertex_->geometric_normal};
                Spectrum pdf_nee = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, light_point, nee_p_cache, scene) *
                    nee_pdf_cache;
                Spectrum dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                // w = make_zero_spectrum();
                // w = 1;
            }
            radiance += current_path_throughput * w * L;
            break;
        }
        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }
        if (!scatter) {
            if (vertex_->material_id == -1) {
                update_medium_id(*vertex_, ray, curr_medium_id);
                continue;
            } 
            never_scatter = false;
            Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            Material mat = scene.materials[vertex_->material_id];
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                            dir_view,
                            *vertex_,
                            scene.texture_pool,
                            bsdf_rnd_param_uv,
                            bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                // BSDF sampling failed. Abort the loop.
                break;
            }
            next_dir = bsdf_sample_->dir_out;
            if (dot(vertex_->geometric_normal, next_dir) * dot(vertex_->geometric_normal, ray.dir) > 0) {
                update_medium_id(*vertex_, ray, curr_medium_id);
            }
            Spectrum f = eval(mat, dir_view, next_dir, *vertex_, scene.texture_pool);
            Real bsdf_pdf = pdf_sample_bsdf(mat, dir_view, next_dir, *vertex_, scene.texture_pool);
            current_path_throughput *= f / bsdf_pdf;
            multi_trans_pdf = make_const_spectrum(1);
            dir_pdf = bsdf_pdf;
        } else {
            never_scatter = false;
            Medium medium = scene.media[curr_medium_id];
            PhaseFunction rho = get_phase_function(medium);
            Vector2 rho_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            auto next_dir_ = sample_phase_function(rho, dir_view, rho_rnd_param_uv);
            if (next_dir_.has_value()) {
                next_dir = next_dir_.value();
                dir_pdf = pdf_sample_phase(rho, dir_view, next_dir);
                multi_trans_pdf = make_const_spectrum(1);
                current_path_throughput *= 
                    eval(rho, dir_view, next_dir) /
                    dir_pdf *
                    get_sigma_s(medium, next_org);
            }
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
    }
    if (x == xx && y == yy) {
        std::cout<< "out: " << radiance << std::endl;
    }
    return radiance;
}
