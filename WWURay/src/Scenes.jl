module Scenes

export HitRecord, Sphere, Scene, TriangleMesh, ray_intersect, create_triangles. Portal
#export has_uvs, has_normals, get_vertex, get_uv, get_normal

using LinearAlgebra

#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..OBJMeshes
using ..Materials



#####################################
###### Generic Scene Data Type ######
#####################################
struct Scene
    background::RGB{Float32}
    objects::Array{Any,1}
    lights::Array{Any,1}
end

""" Structure to store data about an intersection
of a ray with an object (a "hit")."""
mutable struct HitRecord
    t::Float64
    intersection::Vec3
    normal::Vec3
    uv::Union{Vec2,Nothing}
    object
end

# Abstract ray-object intersection function:
# Needs to be implemented for each type of object to be rendered
""" Intersect a ray with an object.
Returns a HitRecord with info about the intersetion, or nothing if
the ray doesn't intersect with the object. """
function ray_intersect(ray::Ray, object) end


##################
##### Sphere #####
##################

# Data type:
struct Sphere
    center::Vec3
    radius::Float64
    material::Material
end


struct Portal
    center::Vec3
    radiusX::Float64
    radiusY::Float64
    material::Material
end

""" Ray-sphere intersection. """
function ray_intersect(ray::Ray, object::Sphere)
    # define constants to use later
    cent = object.center
    d = ray.direction
    p = ray.origin
    r = object.radius
    
    # find the value of t
    a = dot(d,d)
    b = dot(p - cent,d)
    c = dot(p-cent, p-cent)-r^2
    discrim = b*b-a*c
    dx = p[1] - cent[1]
    dy = p[2] - cent[2]
    dz = p[3] - cent[3]

    # not solns if discriminant is < 0
    if discrim < 0
        return nothing
    else
        num = -b-sqrt(discrim)
        if num > 0
            t = num/a
        else
            num = -b+sqrt(discrim)
            if num > 0
                t = num/a
            else
                return nothing
            end
        end
    end
    # intersection points
    x = p[1] + t*d[1]
    y = p[2] + t*d[2]
    z = p[3] + t*d[3]

    # to use for the uv coords, intersection - center
    dist = normalize(Vec3(x,y,z)-cent)

    normal = normalize((Vec3(x,y,z)- cent)/r)
    if isnothing(object.material.texture)
        return HitRecord(t,Vec3(x,y,z), normal, nothing, object)
    else
        theta = atan(dist[1],dist[3])
        phi = asin(-dist[2]) # negative as image was upside down 
        # shift by .5 for middle of the texel
        u = .5 + (theta/ (2*π))
        v = .5 - (phi/π)
        uv = Vec2(u,v)
        return HitRecord(t,Vec3(x,y,z), normal, uv, object)
    end
end


###########################
###### Triangle Mesh ######
###########################

""" Data type: stores the OBJTriangle, a reference to its Mesh
object, and the material it should be rendered with. """
struct Triangle
    geometry::OBJTriangle
    mesh::OBJMesh
    material
end

""" Return an Array of all Triangles belonging to the given mesh, assigning
each one the given material. """
function create_triangles(mesh::OBJMesh, material)
    [Triangle(f, mesh, material) for f in mesh.triangles]
end

""" Some helper functions that make for easier access to triangle data: """
function get_vertex(tri::Triangle, i)
    tri.mesh.positions[tri.geometry.positions[i]]
end

function has_uvs(tri::Triangle)
    length(tri.geometry.uvs) == 3
end

function get_uv(tri::Triangle, i)
    tri.mesh.uvs[tri.geometry.uvs[i]]
end

function has_normals(tri::Triangle)
    length(tri.geometry.normals) == 3
end

function get_normal(tri::Triangle, i)
    tri.mesh.normals[tri.geometry.normals[i]]
end


function ray_intersect(ray::Ray, object::Triangle)



    # vertices
    x = get_vertex(object, 1)
    y = get_vertex(object, 2)
    z = get_vertex(object, 3)
 
    dir = ray.direction
    org = ray.origin
    
    # letters from the text
    a = x[1]-y[1]
    b = x[2]-y[2]
    c = x[3]-y[3]
    g = dir[1]
    h = dir[2]
    i = dir[3]
    j = x[1] - org[1] 
    k = x[2] - org[2]
    l = x[3] - org[3]
    d = x[1]-z[1]
    e = x[2]-z[2]
    f = x[3]-z[3]
    
    # the main matrix
    m = (a * ((e * i) - (h * f))) + (b * ((g * f) - (d * i))) + (c * ((d * h) - (e * g)))
    if m == 0 # no solns
        return nothing
    end
    # cramers rule
    t = -(((f * (a * k - j * b)) + (e * (j * c - a * l)) + (d * (b * l - k * c))) / m)
    γ = ((i * (a * k - j * b)) + (h * (j * c - a * l)) + (g * (b * l - k * c))) / m
    if (γ > 1) || (γ < 0)
        return nothing
    end
    β = (j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g)) / m
    if (β < 0) || (β > (1 -γ))
        return nothing
    end
    # intersection point
    dir *= t
    inter = org + dir
    if !has_normals(object)
        norm = normalize(cross(y-x, z-x))
    else
        # barycentric interpolation for normals
        n1 = get_normal(object, 1)
        n2 = get_normal(object, 2)
        n3 = get_normal(object, 3)
        x = n1[1]+(β*(n2[1]-n1[1]))+(γ*(n3[1]-n1[1]))
        y = n1[2]+(β*(n2[2]-n1[2]))+(γ*(n3[2]-n1[2]))
        z = n1[3]+(β*(n2[3]-n1[3]))+(γ*(n3[3]-n1[3]))
        norm = normalize(Vec3(x,y,z))
    end
    if isnothing(object.material.texture)
        return HitRecord(t, inter, norm, nothing, object)
    else
        if has_uvs(object)
            uv1 = get_uv(object, 1)
            uv2 = get_uv(object, 2)
            uv3 = get_uv(object, 3)

            # barycentric interpolation
            u = uv1[1]+(β*(uv2[1]-uv1[1]))+(γ*(uv3[1]-uv1[1]))
            v = uv1[2]+(β*(uv2[2]-uv1[2]))+(γ*(uv3[2]-uv1[2]))
            uvs = Vec2(u,v)
            return HitRecord(t, inter, norm, uvs, object)
        end
    end
end
end # module Scenes
