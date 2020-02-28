module Materials

using Images
using FileIO

using ..GfxBase

export Material, Texture
export ShadingModel, Flat, Normal
export PhysicalShadingModel, Lambertian, BlinnPhong

export get_diffuse

## Shading model type hierarchy ##
# Types are as follows:
# 
# ShadingModel
#   - Flat - simply the diffuse color of the object
#   - Normal - color-code a pixel according to its surface normal
#   - PhysicalShadingModel
#       - Lambertian - diffuse shading
#       - BlinnPhong - diffuse+specular shading

abstract type ShadingModel end
struct Flat<:ShadingModel end
struct Normal <: ShadingModel end

abstract type PhysicalShadingModel <: ShadingModel end

struct Lambertian <: PhysicalShadingModel end

mutable struct BlinnPhong <: PhysicalShadingModel
    specular_color::RGB{Float32} # color of the highlight
    specular_exp::Float64 # "sharpness" exponent
end

## Texture struct definition ##
mutable struct Texture
    image_data::Array{RGB{Float32},2}
    repeat::Bool
end

""" Texture constructor - loads image data from a file"""
function Texture(image_fn::String, repeat::Bool)
    image_data = load(image_fn)
    Texture(image_data, repeat)
end

## Material struct definition ##
mutable struct Material
    shading_model::ShadingModel
    mirror_coeff::Float64
    texture::Union{Texture,Nothing}
    diffuse_color::Union{RGB{Float32},Nothing}
end

""" Get the diffuse color of a material; if the material is textured,
provide the uv coordinate on the object of that material. """
function get_diffuse(material::Material, uv::Union{Vec2,Nothing})
    if isnothing(material.texture) || isnothing(uv) 
        return material.diffuse_color
    else
        return get_texture_value(material.texture, uv)
    end
end

""" Look up a texture value given the uv coordinates """
function get_texture_value(texture::Texture, uv::Vec2)
    # this equation is from the book
    i = convert(Int, round(uv[1]*size(texture.image_data, 1)-.5))
    j = convert(Int, round(size(texture.image_data,2)-uv[2]*size(texture.image_data,2)-.5))
    i =  i%size(texture.image_data,1)
    j =  j%size(texture.image_data,2)
    
    # make sure we dont index a value that isnt valid
    if i < 1
        i += size(texture.image_data,1)
    end
    if j < 1
        j += size(texture.image_data, 2)
    end
    # get rgb value from uv coord
    rgb = texture.image_data[j, i]
    return rgb
end

end # module Materials
