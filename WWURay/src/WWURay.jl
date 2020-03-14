""" 
Main module for CS480/580 A2 raytracer. Contains core raytracing algrithm,
while referencing several other modules that encapsulate supporting
functionality.
"""

module WWURay

export main

using FileIO
using Images
using StaticArrays
using LinearAlgebra

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
include("Lights.jl")
include("Materials.jl")
include("OBJMeshes.jl")
include("Scenes.jl")
include("Cameras.jl")
include("TestScenes.jl")

using .GfxBase
using .Lights
using .Materials

import .Scenes
import .Scenes.Scene
import .Scenes.HitRecord
import .Cameras
import .TestScenes

# Ray-Scene intersection:
""" Find the closest intersection point among all objects in the scene
along a ray, constraining the search to values of t between tmin and tmax. """
function closest_intersect(objects::Array{Any, 1}, ray::Ray, tmin, tmax)
    closestHit = nothing 
    # various cases for error handling
    for i ∈ objects 
        hit = Scenes.ray_intersect(ray, i)
        if isnothing(hit) && isnothing(closestHit)
            continue
        elseif !isnothing(hit) && isnothing(closestHit)
            if in_bounds(hit.t, tmin, tmax)
                closestHit = hit
            else
                continue
            end
        elseif isnothing(hit) && !isnothing(closestHit)
            if in_bounds(closestHit.t, tmin, tmax)
                closestHit = closestHit
            else
                continue
            end
        elseif hit.t < closestHit.t
            if in_bounds(hit.t, tmin, tmax)
                closestHit = hit
            else
                continue
            end
        end
    end
    return closestHit
end

# helper function to check if a t value is valid
function in_bounds(t, tmin, tmax)
    return (t >= tmin) && (t <= tmax)
end

""" Trace a ray from orig along ray through scene, using Whitted recursive raytracing 
limited to rec_depth recursive calls. """
function traceray(scene::Scene, ray::Ray, tmin, tmax, rec_depth=1)
    # find the closest intersection
    closest_hitrec = closest_intersect(scene.objects, ray, tmin, tmax)
    if isnothing(closest_hitrec)
        return scene.background
    end
    # hitrec vlaues
    object = closest_hitrec.object
    point = closest_hitrec.intersection
    normal = closest_hitrec.normal
    material = object.material
    shader = material.shading_model

    # determine the color of the pixel
    local_color = determine_color(shader, object.material, ray, closest_hitrec, scene)
    mirror_co = object.material.mirror_coeff
    # ray direction as per the lecture
    ray_dir = 2* normal*dot(normal, -ray.direction)+ray.direction 
    mirror_ray = Ray(point, normalize(ray_dir))

    # recurse until rec depth is hit or if there is no mirror coef
    if (mirror_co == 0) || (rec_depth <= 0)
        return local_color
    end
    reflected_color = traceray(scene, mirror_ray, 1e-8, Inf, rec_depth-1)
    return (local_color*(1-mirror_co))+(reflected_color*mirror_co)
        
end

""" Determine the color of interesction point described by hitrec 
Flat shading - just color the pixel the material's diffuse color """
function determine_color(shader::Flat, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    get_diffuse(material, hitrec.uv)
end
""" Normal shading - color-code pixels according to their normals """
function determine_color(shader::Normal, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    normal_color = normalize(hitrec.normal) / 2 .+ 0.5
    RGB{Float32}(normal_color...)
end


""" Determine the color of a physical (Lambertian, BlinnPhong, etc.) surface """
function determine_color(shader::PhysicalShadingModel, material::Material, ray::Ray, hitrec::HitRecord, scene::Scene)
    # start with black color, increment if point is not shadowed
    lightAmount = RGB(0.0,0.0,0.0)
    for i ∈  scene.lights
        # call light source depending on type
        if !is_shadowed(scene, i, hitrec.intersection)
            lightAmount += shade_light(shader, material, ray, hitrec, i, scene) 
        end
    end
    return lightAmount

end

""" shade_light(shader, material, ray, hitrec, light, scene)
Determine the color contribution of the given light along the given ray.
Color depends on the material, the shading model (shader), properties of the intersection 
given in hitrec, """
# just diffuse shading
function shade_light(shader::Lambertian, material::Material, ray::Ray, hitrec, light, scene)
    l = get_diffuse(material, hitrec.uv) * light.intensity*max(0, dot(hitrec.normal, normalize(Lights.light_direction(light, hitrec.intersection))))
end

""" Blinn-Phong surface shading """
function shade_light(shader::BlinnPhong, material::Material, ray::Ray, hitrec, light, scene)
    # diffuse coefficent
    diffuse = get_diffuse(material, hitrec.uv) * light.intensity*max(0, dot(hitrec.normal, normalize(Lights.light_direction(light, hitrec.intersection))))
    # half vector
    h = normalize((-ray.direction)+Lights.light_direction(light, hitrec.intersection))
    # specular reflection
    spec = shader.specular_color*light.intensity*(max(0, dot(hitrec.normal, h))^(shader.specular_exp))
    return diffuse+spec
end


""" Determine whether point is in shadow wrt light """
function is_shadowed(scene, light::DirectionalLight, point::Vec3)
    if closest_intersect(scene.objects, Ray(point, Lights.light_direction(light,point)), 1e-8, Inf) != nothing
        return true
    else
        return false
    end
end

function is_shadowed(scene, light::PointLight, point::Vec3)
    if closest_intersect(scene.objects, Ray(point, Lights.light_direction(light, point)), 1e-8, 1) != nothing
        return true
    else
        return false
    end
end

# Main loop:
function main(scene, height, width, out)
    if scene == 1
	    scene = TestScenes.portalScene1(height, width)[1]
        camera = TestScenes.portalScene1(height, width)[2]
        getImages(scene, camera, height, width, out, 1, 0)
    elseif scene == 2
	    scene = TestScenes.portalScene2(height, width)[1]
        camera = TestScenes.portalScene2(height, width)[2]
        getImages(scene, camera, height, width, out, 1, -1)
    elseif scene == 3
	    scene = TestScenes.portalScene3(height, width)[1]
        camera = TestScenes.portalScene3(height, width)[2]
        getImages(scene, camera, height, width, out, 1, 0)   
   end
end


function getImages(scene, camera, height, width, out, amountx, amountz) 

    RANGLE = 45         # Rotate angle size
    NUMS = 60           # Number of images
    
    for idx = 1:NUMS
    
        # change camera
        T = [
            1.0 0.0 0.0 inv(tan(RANGLE/NUMS))/30 * amountx;
            0.0 1.0 0.0 0.0;
            0.0 0.0 1.0 inv(tan(RANGLE/NUMS))/30 * amountz;
            0.0 0.0 0.0 1.0
        ] 
        
        R = TestScenes.rotate(RANGLE/NUMS, 0, 1, 0);
        newcam = T * Vec4(camera.eye[1], camera.eye[2], camera.eye[3], 1)
        camera.eye = Vec3(newcam[1], newcam[2], newcam[3])


        #get objects
        objs = scene.objects
        # Create a blank canvas to store the image:
        canvas = zeros(RGB{Float32}, height, width)

        # for each pixel iterate over each object
        # color the canavs accordingly 
        for i = 1:height
            for j = 1:width
                ray:: Ray = Cameras.pixel_to_ray(camera, i, j)
                color = traceray(scene, ray, 1, Inf, 8)
                canvas[i,j] = color
            end
        end

        # clamp canvas to valid range:
        clamp01!(canvas)
        
        outfolder = string("video/", out)
        outfname = string(idx, ".png")
        save(File(format"PNG", string(outfolder, outfname)), colorview(RGB, canvas))

    end # end of idx loop    
end


end # module WWURay