#include "scene.h"
#include "intersect.h"
#include "montecarlo.h"
#include "animation.h"

#include <thread>

using std::thread;

// modify the following line to disable/enable parallel execution of the pathtracer
bool parallel_pathtrace = true;

image3f pathtrace(Scene* scene, bool multithread);
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose);
vec3f pathtrace_ray(Scene* scene, ray3f ray, Rng* rng, int depth);

//  texture value
vec3f lookup_scaled_texture(vec3f value, image3f* texture, vec2f uv, intersection3f intersection,
                            float distance = 0.0, bool tile = false, bool bilinear_filter = false, bool mipmap_filter = false) {

    if(not texture) return value;

    // Declare common variables
    float u, v;

    // if tiling, then tile
    if (tile) {
        // Tile the texture if texcoords exceed [0, 1)
        u = uv.x - floor(uv.x);
        v = uv.y - floor(uv.y);
    } else {
        // Give Clamp coordinates
        u = clamp(uv.x, 0.0f, 1.0f);
        v = clamp(uv.y, 0.0f, 1.0f);
    }

    // Return the simple value from one pixel if not bilinear or trilinear filtering (mipmap)
    if (!bilinear_filter && !mipmap_filter) {
        return value * texture->at(u*(texture->width()-1), v*(texture->height()-1));
    }

    // Arrange the textures into an indexable array
    vector<image3f*> mip_maps;
    mip_maps.push_back(intersection.mat->mipmap_texture_1);
    mip_maps.push_back(intersection.mat->mipmap_texture_2);
    mip_maps.push_back(intersection.mat->mipmap_texture_3);
    mip_maps.push_back(intersection.mat->mipmap_texture_4);
    mip_maps.push_back(intersection.mat->mipmap_texture_5);
    mip_maps.push_back(intersection.mat->mipmap_texture_6);
    mip_maps.push_back(intersection.mat->mipmap_texture_7);
    mip_maps.push_back(intersection.mat->mipmap_texture_8);
    mip_maps.push_back(intersection.mat->mipmap_texture_9);
    mip_maps.push_back(intersection.mat->mipmap_texture_10);


    // Compute interpolation factors
    int texture_count = mip_maps.size() - 1;                        // number of different resolution texture image
    int i = (int)(u * texture->width());                                // slides on texturing
    int j = (int)(v * texture->height());

    // Declare filtered color vector
    vec3f interpolated_color;

    // if bilinear, then bilinear filter
    if (!mipmap_filter) {

        float s = texture->width() * u - floor(texture->width() * u);
        float t = texture->height() * v - floor(texture->height() * v);

        // Retreive the pixels that we will interpolate with
        vec3f c_0_0 = texture->at(i-1, j-1);
        vec3f c_0_1 = texture->at(i, j-1);
        vec3f c_1_0 = texture->at(i-1, j);
        vec3f c_1_1 = texture->at(i, j);

        // Find the interpolated color vector
        interpolated_color = (
                    (1 - s) * (1 - t) * c_0_0 +
                    (s) * (1 - t) * c_0_1 +
                    (1 - s) * (t) * c_1_0 +
                    s * t * c_1_1
        );
        return value * interpolated_color;
    }

    // if mip map, then mip map

    // Calculate Interpolation Factors
    float distance_heuristic = clamp( ((distance - 1) / 5), 0.0, 0.8);
    float s = texture->width() * u - (int)(texture->width() * u);
    float t = texture->height() * v - (int)(texture->height() * v);
    float w = texture_count * distance_heuristic - (int)(texture_count * distance_heuristic);
    int k = (int)(distance_heuristic * texture_count);

    // Get the pixels
    vec3f c_0_0_0 = mip_maps[k]->at(i, j);
    vec3f c_1_0_0 = mip_maps[k]->at(i+1, j);
    vec3f c_0_1_0 = mip_maps[k]->at(i, j+1);
    vec3f c_0_0_1 = mip_maps[k+1]->at(i, j);
    vec3f c_1_0_1 = mip_maps[k+1]->at(i+1, j);
    vec3f c_0_1_1 = mip_maps[k+1]->at(i, j+1);
    vec3f c_1_1_0 = mip_maps[k]->at(i+1, j+1);
    vec3f c_1_1_1 = mip_maps[k+1]->at(i+1, j+1);

    // INterpolate the new color using trilinear interpolation (shirley 713)
    interpolated_color = (
                (1-s) *   (1 - t) *   (1 - w) *   c_0_0_0 +
                (1-s) *   (1 - t) *   w     *   c_0_0_1 +
                (1-s) *   t     *   (1 - w) *   c_0_1_0 +
                s     *   (1 - t) *   (1 - w) *   c_1_0_0 +
                s     *   (1 - t) *   w     *   c_1_0_1 +
                (1-s) *   t     *   w     *   c_0_1_1 +
                s     *   t     *   (1 - w) *   c_1_1_0 +
                s     *   t     *   w     *   c_1_1_1
    );
    return value * interpolated_color;


}

// evaluate the environment map
vec3f eval_env(vec3f ke, image3f* ke_txt, vec3f dir, intersection3f intersection) {
    auto u = atan2(dir.x, dir.z)/(2 * pif);
    auto v = 1 - acos(dir.y)/(pif);
    auto env = lookup_scaled_texture(ke, ke_txt, vec2f(u, v), intersection ,0.0, true, false, false);
    return env;
}


// compute the brdf
vec3f eval_brdf(vec3f kd, vec3f ks, float n, vec3f v, vec3f l, vec3f norm, bool microfacet) {
    if (not microfacet) {
        auto h = normalize(v+l);
        return kd/pif + ks*(n+8)/(8*pif) * pow(max(0.0f,dot(norm,h)),n);
    } else {
        auto h = normalize(v+l);
        auto d = (n+2)/(2 * pif) * pow(max(0.0, dot(h, norm)), n);
        auto f = ks + (one3f - ks) * pow((1 - dot(h, l)), 5);
        auto g = min(min(1.0f, 2*dot(h,norm)*dot(v,norm)/dot(v,h)), 2*dot(h,norm)*dot(l,norm)/dot(v,h));
        auto p = d*g*f/(4*dot(l, norm)*dot(v, norm));
        return p;
    }
}


// compute the color corresponing to a ray by pathtrace
vec3f pathtrace_ray(Scene* scene, ray3f ray, Rng* rng, int depth) {

    // get scene intersection
    auto intersection = intersect(scene,ray);
    bool mipmap_filter;

    // if not hit, return background (looking up the texture by converting the ray direction to latlong around y)
    if(not intersection.hit) {
        return eval_env(scene->background, scene->background_txt, ray.d, intersection);
    } else {
        mipmap_filter = intersection.mat->mipmap_filter;
    }

    // setup variables for shorter code
    auto pos = intersection.pos;
    auto norm = intersection.norm;
    auto v = -ray.d;

    vec3f kd_1, kd_2, kd_3;

    bool tiling = intersection.mat->tiling;
    bool bilinear_filter = intersection.mat->bilinear_filter;

    // compute material values by looking up textures
    float distance = dist(intersection.pos, scene->camera->frame.o);
    auto kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, intersection.texcoord, intersection, distance, tiling, bilinear_filter, mipmap_filter);
    auto ke = lookup_scaled_texture(intersection.mat->ke, intersection.mat->ke_txt, intersection.texcoord, intersection, distance, tiling, bilinear_filter, false);
    auto ks = lookup_scaled_texture(intersection.mat->ks, intersection.mat->ks_txt, intersection.texcoord, intersection, distance, tiling, bilinear_filter, false);
    auto n = intersection.mat->n;
    auto mf = intersection.mat->microfacet;

    // accumulate color starting with ambient
    auto c = scene->ambient * kd;

    // add emission if on the first bounce
    if(depth == 0 and dot(v,norm) > 0) c += ke;
    // foreach point light
     for(auto light : scene->lights) {
         // compute light response
         auto cl = light->intensity / (lengthSqr(light->frame.o - pos));
         // compute light direction
         auto l = normalize(light->frame.o - pos);
         // compute the material response (brdf*cos)
         auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, mf);
         // multiply brdf and light
         auto shade = cl * brdfcos;
         // check for shadows and accumulate if needed
         if(shade != zero3f){
             // if shadows are enabled
             if(scene->path_shadows) {
                 // perform a shadow check and accumulate
                 if(not intersect_shadow(scene,ray3f::make_segment(pos,light->frame.o))) c += shade;
             } else {
                 // else just accumulate
                 c += shade;
             }
         }
     }

     // foreach surface
     for(Surface* surface : scene->surfaces){
         // skip if no emission from surface
         if (surface->mat->ke != zero3f){
             // todo: pick a point on the surface, grabbing normal, area, and texcoord
             vec2f ruv;
             float sa;
             vec3f snorm;
             vec3f lpos;
             // check if quad
             if (surface->isquad){
                 // generate a 2d random number
                 ruv = rng->next_vec2f();
                 // compute light position, normal, area
                 lpos = vec3f(surface->radius * (2 * ruv.x - 1), surface->radius * (2 * ruv.y - 1), 0);
                 lpos = transform_point_from_local(surface->frame, lpos);
                 sa = surface->radius * surface->radius;
                 snorm = transform_vector_from_local(surface->frame, vec3f(0, 0, 1));
                 // set tex coords as random value got before

             }
             // else if sphere
             else{
                 // generate a 2d random number
                 ruv = rng->next_vec2f();

                 // compute light position, normal, area
                 lpos = sample_direction_spherical_uniform(ruv) * surface->radius / length(sample_direction_spherical_uniform(ruv));
                 snorm = normalize(lpos);
                 snorm = transform_vector_from_local(surface->frame, snorm);
                 lpos = transform_point_from_local(surface->frame, lpos);
                 sa = 4 * PI * pow(surface->radius, 2);
                 // set tex coords as random value got before

             }

             // get light emission from material and texture
             // compute light direction
             vec3f l = normalize(lpos - pos);
             // compute light response (ke * area * cos_of_light / dist^2)
             auto cl = (lookup_scaled_texture(surface->mat->ke, surface->mat->ke_txt, ruv, intersection, false) * sa * 4 * max(-dot(snorm, l), 0.0)/distSqr(pos, lpos));
             // compute the material response (brdf*cos)
             auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, mf);
             // multiply brdf and light
             auto shade = brdfcos * cl;
             // check for shadows and accumulate if needed
             if(shade != zero3f){
                 // if shadows are enabled
                 if(scene->path_shadows) {
                     // perform a shadow check and accumulate
                     if(not intersect_shadow(scene,ray3f::make_segment(pos,lpos))) c += shade;
                 }
                 // else just accumulate
                 else{
                     c += shade;
                 }
             }
         }
     }

     // todo: sample the brdf for environment illumination if the environment is there
     // if scene->background is not zero3f
     if (scene->background != zero3f){
         // pick direction and pdf
         //vec2f ruv = vec2f((rand()/(float)(RAND_MAX + 1)), (rand()/(float)(RAND_MAX + 1)));
         vec2f ruv = rng->next_vec2f();
         auto lpdf = sample_brdf(kd, ks, n, v, norm, ruv, rng->next_float());
         vec3f l = lpdf.first;
         float pdf = lpdf.second;
         // compute the material response (brdf*cos)
         auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, mf);
         // todo: accumulate response scaled by brdf*cos/pdf
         auto cl = eval_env(scene->background, scene->background_txt, l, intersection);
         auto shade = cl * brdfcos/pdf;
         // if material response not zero3f
         if(shade != zero3f){
             // if shadows are enabled
             if(scene->path_shadows){
                 // perform a shadow check and accumulate
                 if(not intersect_shadow(scene,ray3f(pos, l)))
                     c += shade;
             }
             // else just accumulate
             else{
                 c += shade;
             }
         }
     }

     // todo: sample the brdf for indirect illumination
     // if kd and ks are not zero3f and haven't reach max_depth
     if ((kd != zero3f || ks !=zero3f) && depth < scene->path_max_depth){
         // pick direction and pdf
         vec2f ruv = rng->next_vec2f();
         auto lpdf = sample_brdf(kd, ks, n, v, norm, ruv, rng->next_float());
         vec3f l = lpdf.first;
         float pdf = lpdf.second;
         // compute the material response (brdf*cos)
         auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, mf);
         // accumulate recersively scaled by brdf*cos/pdf
         c += brdfcos/pdf * pathtrace_ray(scene, ray3f(pos, l),rng, depth+1);
     }


     // if the material has reflections
        if( !(intersection.mat->kr == zero3f)) {
            // create the reflection ray
            auto reflection_ray = ray3f(intersection.pos,reflect(ray.d,intersection.norm));
                   // if there is a blur
                   if(intersection.mat->blurred_reflection_size != 0 && intersection.mat->num_blurred_samples != 0) {
                       vec3f refl_ray_i = zero3f;
                       vec3f u = vec3f((float)(-reflection_ray.d.y), reflection_ray.d.x, 0);
                       vec3f v = vec3f(0, (float)(-reflection_ray.d.z), reflection_ray.d.y);

                       for (int i = 0; i < intersection.mat->num_blurred_samples; i++) {
                           vec2f random_ray = rng->next_vec2f();
                           vec3f new_vector = reflection_ray.d + ((float)(0.5) - random_ray.x) * intersection.mat->blurred_reflection_size * u +
                                                      ((float)0.5 - random_ray.y) * intersection.mat->blurred_reflection_size * v;

                           new_vector = normalize(new_vector);
                           ray3f new_ray = ray3f(reflection_ray.e, new_vector );
                           refl_ray_i += pathtrace_ray(scene, new_ray, rng, depth + 1);
                       }
                        c += (intersection.mat->kr / (float)intersection.mat->num_blurred_samples) * refl_ray_i;

                   } else {
                       // get all reflected light
                       c += pathtrace_ray(scene,reflection_ray,rng,depth+1);
                   }
        }

     // return the accumulated color
     return c;

}

// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "04_pathtrace", "raytrace a scene",
            {  {"resolution",     "r", "image resolution", typeid(int),    true,  jsonvalue() } },
            {  {"scene_filename", "",  "scene filename",   typeid(string), false, jsonvalue("scene.json") },
               {"image_filename", "",  "image filename",   typeid(string), true,  jsonvalue("") } }
        });

    auto scene_filename = args.object_element("scene_filename").as_string();
    Scene* scene = nullptr;
    if(scene_filename.length() > 9 and scene_filename.substr(0,9) == "testscene") {
        int scene_type = atoi(scene_filename.substr(9).c_str());
        scene = create_test_scene(scene_type);
        scene_filename = scene_filename + ".json";
    } else {
        scene = load_json_scene(scene_filename);
    }
    error_if_not(scene, "scene is nullptr");

    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";

    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }

    // NOTE: acceleration structure does not support animations
    message("reseting animation...\n");
    animate_reset(scene);

    message("accelerating...\n");
    accelerate(scene);

    message("rendering %s...\n", scene_filename.c_str());
    auto image = pathtrace(scene, parallel_pathtrace);

    message("saving %s...\n", image_filename.c_str());
    write_png(image_filename, image, true);

    delete scene;
    message("done\n");
}


/////////////////////////////////////////////////////////////////////
// Rendering Code


// pathtrace an image
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose) {
    if(verbose) message("\n  rendering started        ");
    // foreach pixel
    for(auto j = offset_row; j < scene->image_height; j += skip_row ) {
        if(verbose) message("\r  rendering %03d/%03d        ", j, scene->image_height);
        for(auto i = 0; i < scene->image_width; i ++) {
            // init accumulated color
            image->at(i,j) = zero3f;
            // grab proper random number generator
            auto rng = &rngs->at(i, j);
            // foreach sample
            for(auto jj : range(scene->image_samples)) {
                for(auto ii : range(scene->image_samples)) {
                    // compute ray-camera parameters (u,v) for the pixel and the sample
                    auto u = (i + (ii + rng->next_float())/scene->image_samples) /
                        scene->image_width;
                    auto v = (j + (jj + rng->next_float())/scene->image_samples) /
                        scene->image_height;
                    // compute camera ray
                    auto ray = transform_ray(scene->camera->frame,
                        ray3f(zero3f,normalize(vec3f((u-0.5f)*scene->camera->width,
                                                     (v-0.5f)*scene->camera->height,-1))));
                    // set pixel to the color raytraced with the ray
                    image->at(i,j) += pathtrace_ray(scene,ray,rng,0);
                }
            }
            // scale by the number of samples
            image->at(i,j) /= (scene->image_samples*scene->image_samples);
        }
    }
    if(verbose) message("\r  rendering done        \n");

}

// pathtrace an image with multithreading if necessary
image3f pathtrace(Scene* scene, bool multithread) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);

    // create a random number generator for each pixel
    auto rngs = RngImage(scene->image_width, scene->image_height);

    // if multitreaded
    if(multithread) {
        // get pointers
        auto image_ptr = &image;
        auto rngs_ptr = &rngs;
        // allocate threads and pathtrace in blocks
        auto threads = vector<thread>();
        auto nthreads = thread::hardware_concurrency();
        for(auto tid : range(nthreads)) threads.push_back(thread([=](){
            return pathtrace(scene,image_ptr,rngs_ptr,tid,nthreads,tid==0);}));
        for(auto& thread : threads) thread.join();
    } else {
        // pathtrace all rows
        pathtrace(scene, &image, &rngs, 0, 1, true);
    }

    // done
    return image;
}


