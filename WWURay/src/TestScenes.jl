module TestScenes

#push!(LOAD_PATH, pwd())
using ..GfxBase
using ..Scenes
using ..Materials
using ..Lights
using ..OBJMeshes
using ..Cameras
using LinearAlgebra

# helpful things:
make_diffuse(color) = Material(Lambertian(), 0.0, nothing, color)
black = RGB{Float32}(0,0,0)
red = RGB{Float32}(1,0,0)
green = RGB{Float32}(0,1,0)
blue = RGB{Float32}(0,0,1)
white = RGB{Float32}(1,1,1)
purple = RGB{Float32}(1,0,1)

function camera_1(img_height, img_width)
    CanonicalCamera(img_height, img_width)
end



function camera_portal1(img_height, img_width)
    eye = Vec3(-3.5, 0.5, 5)
    view = Vec3(-2.5, 0.5, 0) - eye
    up = Vec3(0, 1, 0)
    focal = 1.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end




function camera_portal2(img_height, img_width)
    eye = Vec3(-6, 0.5, 8)
    view = Vec3(0, 0, 2) - eye
    up = Vec3(0, 1, 0)
    focal = 1.0
    Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end






function camera_portal3(img_height, img_width)
	eye = Vec3(-0,0.5,5)
	view = Vec3(-0.2,0.45,3) - eye
	up = Vec3(0,1,0)
	focal = 1.0
	Cameras.PerspectiveCamera(eye, view, up, focal, img_height, img_width)
end






cameras = [camera_1, camera_portal1, camera_portal2, camera_portal3]

function get_camera(i, img_height, img_width)
    cameras[i](img_height, img_width)
end


function get_scene(i)
    scenes[i]()
end

function scene_1()
    bg = RGB{Float32}(0.95, 0.95, 0.95)
    objs = [Sphere(Vec3(0, 0, -5), 1, Material(Flat(), 0.0, nothing, RGB{Float32}(0.73,0,0.17)))]
    lights = [PointLight(0.8, Vec3(0,0,0))]
    Scene(bg, objs, lights)
end

function scene_2()
    bg = black
    objs = [
            Sphere(Vec3( 2, 0, -8), 1, Material(Lambertian(), 0.0, nothing, white)),
            Sphere(Vec3(-2, 0, -8), 2, Material(Lambertian(), 0.0, nothing, blue))
           ]

    lights = [ DirectionalLight(1.0, Vec3(1, 0.5, -0.1)) ]
    Scene(bg, objs, lights)
end

function scene_3()
    bg = black
    mat = Material(Lambertian(), 0.0, nothing, white)
    objs = [
            Sphere(Vec3( -2, 1, -8), 1, mat),
            Sphere(Vec3(  2, 1, -8), 1, mat)
           ]

    lights = [ PointLight(1.0, Vec3(0, 5, -8.5)) ]

    Scene(bg, objs, lights)
end

function scene_4()
    bg = black
    mat1 = Material(BlinnPhong(white, 10), 0.0, nothing, white)
    mat2 = Material(BlinnPhong(white, 10), 0.0, nothing, blue)
    mat3 = Material(BlinnPhong(red, 100), 0.0, nothing, blue)
    objs = [
            Sphere(Vec3( -2, -1, -8), 1, mat1),
            Sphere(Vec3( -1, 1, -8), 1, mat2),
            Sphere(Vec3(  0, -1, -8), 1, mat3),
            Sphere(Vec3(  1, 1, -8), 1, mat2),
            Sphere(Vec3(  2, -1, -8), 1, mat1),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.0, nothing, white))
           ]

    lights = [ PointLight(0.8, Vec3(0, 4, -8)),
               PointLight(0.2, Vec3(0, 0, 0)) ]

    Scene(bg, objs, lights)
end

function scene_5()
    bg = black

    mat = Material(Lambertian(), 0.0, nothing, white)

    objs = [
            Sphere(Vec3( -1, 0, -6), 0.5, mat),
            Sphere(Vec3(  1, 0, -5), 0.5, Material(Lambertian(), 0.4, nothing, white)),
            Sphere(Vec3( -1, 0, -4), 0.5, mat),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, white)) # ground
           ]

    lights = [ DirectionalLight(0.6, Vec3(1, 1, 0)),
               PointLight(0.4, Vec3(0, 0, 0)) ]

    Scene(bg, objs, lights)
end

function scene_6()
    bg = black

    r = Material(BlinnPhong(white, 10), 0.0, nothing, red)
    g = Material(BlinnPhong(white, 10), 0.0, nothing, green)
    b = Material(BlinnPhong(white, 10), 0.0, nothing, blue)
    refl = Material(Lambertian(), 0.6, nothing, white)

    objs = [
            #Sphere(Vec3(-10, 0, -1), 9.2, refl),
            Sphere(Vec3(-1,  -1.1, -3), 0.5, r),
            Sphere(Vec3( -0.5,  -1.0, -4), 0.5, g),
            Sphere(Vec3( 0,  -0.9, -5), 0.5, b),
            Sphere(Vec3( 5,  -1, -4), 4, refl),
            #Sphere(Vec3( 10,  0.1 , -1), 9.2, refl),
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, white)) # floor
           ]

    lights = [ PointLight(0.6, Vec3(1, 10, -4)),
               PointLight(0.4, Vec3(0, 0, 0)) ]

    Scene(bg, objs, lights)
end



""" Take the OBJMesh mesh and return an array of Triangles from the mesh
with the given material, after scaling the mesh positions by scale and moving
them by translation """
function mesh_helper(mesh, material, scale=1.0, translation=Vec4(0,0,0,0), rotAx=[Vec3(1,0,0)], angle=[0])
    rot = rotate(angle[1], rotAx[1][1], rotAx[1][2], rotAx[1][3])
    for i in 2:length(rotAx)
        rot = rotate(angle[i], rotAx[i][1], rotAx[i][2], rotAx[i][3]) * rot
    end
    for i in 1:length(mesh.positions)
        j = rot*Vec4(mesh.positions[i][1], mesh.positions[i][2], mesh.positions[i][3],1)  * scale + translation
        mesh.positions[i] = Vec3(j[1], j[2], j[3])
    end

    create_triangles(mesh, material)
end

function rotate(a, x, y, z)
    m = Matrix{Float64}(I,4,4)
	ct = cos(a * pi / 180.0);
    st = sin(a * pi / 180.0);
    mag = sqrt(x^2 + y^2 + z^2);
    x = x / mag;
    y = y / mag;
    z = z / mag;

    m[1] = ct + x*x*(1-ct);
    m[2] = x*y*(1-ct) - z*st;
    m[3] = x*z*(1-ct) + y*st;

    m[5] = y*x*(1-ct)+z*st;
    m[6] = ct + y*y*(1-ct);
    m[7] = y*z*(1-ct)-x*st;

    m[9] = z*x*(1-ct)-y*st;
    m[10] = z*y*(1-ct)+x*st;
    m[11] = ct + z*z*(1-ct);
    return m
end





function scene_7()
    bg = black
    objs = []

    # add a bunny:
    bunny_mat = Material(Lambertian(), 0.0, nothing, RGB{Float32}(0.6, 0.5, 0.5))
    bunny = read_obj("data/bunny.obj")
    append!(objs, mesh_helper(bunny, bunny_mat, 1.0, Vec3(0.2, 0, -5)))

    # add a cube
    cube_mat = Material(Lambertian(), 0.6, nothing, white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 10.0, Vec3(-11.2, 0, 0)))

    lights = [ PointLight(0.5, Vec3(1,2,-5)),
               DirectionalLight(0.3, Vec3(0,0,1)),
               DirectionalLight(0.3, Vec3(0,1,1)),
               DirectionalLight(0.3, Vec3(1,1,1)),
               DirectionalLight(0.3, Vec3(0,1,0)) ]

    Scene(bg, objs, lights)

end

scene_8 = scene_7

function scene_9()
    bg = black

    objs = []

    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.2, nothing, RGB{Float32}(0.8, 0.8, 1.0))))
    sphere_material = Material(Lambertian(), 0.0, Texture("data/earth.png", false), nothing)
    push!(objs, Sphere(Vec3(-1.25, 0, -6), 1, sphere_material))

    sphere_m = sphere_mesh(32, 16)
    scale = 1.0
    translation = Vec3(1.25, 0, -6)
    for i in 1:length(sphere_m.positions)
        sphere_m.positions[i] = sphere_m.positions[i] * scale + translation
    end
    append!(objs, create_triangles(sphere_m, sphere_material))

    cube_mat = Material(Lambertian(), 0.0, Texture("data/1.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.5, Vec3(-1, -1, -3)))

    lights = [ DirectionalLight(0.4, Vec3(0,1,0)),
               DirectionalLight(0.8, Vec3(0.4,0.4,1)) ]

    Scene(bg, objs, lights)

end


#artifact for Ryan Wells
function artifact_wellsr8(img_height, img_width)


    bg = RGB{Float32}(0.2,0.2,0.2)
    objs = []

    num = 0.0
    darker = 0.0

    # add bunnies:
    for i = 1:10
        bunny_mat = Material(Lambertian(), 0.0, nothing, RGB{Float32}(1, 0.4 + darker, 0.4 + darker))
        bunny = read_obj("data/bunny.obj")
									      #left/right      up/down      front/back
        append!(objs, mesh_helper(bunny, bunny_mat, 1.0 - (darker / 1.5), Vec3(-1.2 + (num * 7.5), 0 + darker, -5 - (num * 7))))   
        num += 0.1
		darker -= 0.05
    end

    # add a cube
    cube_mat = Material(Lambertian(), 0.40, nothing, black)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 10.0, Vec3(0, -12, 0)))

    # add a cube
    cube_mat = Material(Lambertian(), 0.40, nothing, RGB{Float32}(0.1, 0.01, 0.01))
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 15.0, Vec3(0, 0, -32)))


    lights = [ PointLight(0.5, Vec3(1,2,-5)),
               DirectionalLight(0.3, Vec3(0,0,1)),
               DirectionalLight(0.3, Vec3(0,1,1)),
               DirectionalLight(0.3, Vec3(1,1,1)),
               DirectionalLight(0.3, Vec3(0,1,0)) ]

	lights = [PointLight(0.8, Vec3(0,0,0)),
	PointLight(0.8, Vec3(0,5,0))]

    return (Scene(bg, objs, lights), camera_1(img_height, img_width))

end



function artifact_shinm(img_height, img_width)
    bg = black

    objs = [
            Sphere(Vec3( 0, 2.5, -4), 1.5, Material(Lambertian(), 0.6, nothing, white)),               # reflection
        
            Sphere(Vec3(-1.5, -0.5, -2), 0.4, Material(BlinnPhong(white, 10), 0.0, nothing, red)),
            Sphere(Vec3(-0.5, -0.5, -5), 0.3,  Material(BlinnPhong(white, 10), 0.0, nothing, blue)),
            Sphere(Vec3(0,  -0.5, -5), 0.3, Material(BlinnPhong(white, 10), 0.0, nothing, green)),
            Sphere(Vec3(1, -0.5, -3), 0.4, Material(Lambertian(), 0.0, Texture("data/ball.png", false), nothing)),
            
            Sphere(Vec3(  0, -5001, 0), 5000, Material(Lambertian(), 0.5, nothing, purple)) # floor
           ]

    lights = [ PointLight(0.4, Vec3(-1, 1, 0)), 
               DirectionalLight(1.0, Vec3(1, 0.5, -0.1))]

   return (Scene(bg, objs, lights), camera_2(img_height, img_width))
end


function portalScene1(img_height, img_width)

    bg = black

    objs = []


	#make the ground
    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.02, nothing, RGB{Float32}(0.8, 0.8, 1.0))))

	#central cube
    cube_mat = Material(Lambertian(), 0.0, Texture("data/wall.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 1.0, Vec4(0, 0, 0,1)))
    
    #bunny location
    # add a bunny:
    bunny_mat = Material(Lambertian(), 0.02, nothing, RGB{Float32}(0.6, 0.5, 0.5))
    bunny = read_obj("data/bunny.obj")
    append!(objs, mesh_helper(bunny, bunny_mat, 1.0, Vec4(5, 0, 0,1), [Vec3(0,1,0)], [90]))
    
    
    #background behind bunny
    cube_mat = Material(Lambertian(), 0.0, Texture("data/warningAll.jpg", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 2.0, Vec4(7.5, 0.5, 0,1), [Vec3(0,1,0)], [90]))


    #front portal background
    portal_mat = Material(Lambertian(), 0.0,  nothing, RGB{Float32}(0.99215686274, 0.4, 0.0))
    append!(objs, mesh_helper(portal_background_mesh(2,1), portal_mat, 0.49, Vec4(0,0, 0.52,1), [Vec3(0,0,1)], [90]))

    #front portal
    portal_mat = Material(Lambertian(), 1.0,  nothing, RGB{Float32}(0.0, 0.0, 0.0))
    append!(objs, mesh_helper(portal_mesh(2,1), portal_mat, 0.47, Vec4(0, 0.015, 0.545,1), [Vec3(0,0,1)], [90]))


    #right portal background
    portal_mat = Material(Lambertian(), 0.0,  nothing, RGB{Float32}(0.0, 0.47058823529, 1))
    append!(objs, mesh_helper(portal_background_mesh(2,1), portal_mat, 0.49, Vec4(0.52, 0, 0,1), [Vec3(0,1,0),Vec3(1,0,0)], [-90,90]))

    #right portal
    portal_mat = Material(Lambertian(), 1.0,  nothing, RGB{Float32}(0.0, 0.0, 0.0))
    append!(objs, mesh_helper(portal_mesh(2,1), portal_mat, 0.47, Vec4(0.545, 0.015, 0,1), [Vec3(0,1,0),Vec3(1,0,0)], [-90,90]))


    
    
    lights = [ DirectionalLight(0.8, Vec3(-1, 0.4, 0.8)),
               PointLight(0.9, Vec3(2, 8, 0)),
	       PointLight(0.9, Vec3(3, 3, 3)) ]


   return (Scene(bg, objs, lights), camera_portal1(img_height, img_width))


end




function portalScene2(img_height, img_width)

    bg = black

    objs = []

	#make the ground
    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.03, nothing, RGB{Float32}(0.8, 0.8, 1.0))))

	#left cube
    cube_mat = Material(Lambertian(), 0.0, Texture("data/wall.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 1.0, Vec4(0, 0, 0,1)))

	#right cube
    cube_mat = Material(Lambertian(), 0.0, Texture("data/wall.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 1.0, Vec4(2.005, 0, 2.005,1)))


    #left portal background
    portal_mat = Material(Lambertian(), 0.0,  nothing, RGB{Float32}(0.99215686274, 0.4, 0.0))
	append!(objs, mesh_helper(portal_background_mesh(2,1), portal_mat, 0.49, Vec4(0, 0, 0.515,1), [Vec3(0,0,1)], [90]))

    #left portal
    portal_mat = Material(Lambertian(), 1.0,  nothing, RGB{Float32}(0.0, 0.0, 0.0))
    append!(objs, mesh_helper(portal_mesh(2,1), portal_mat, 0.47, Vec4(0, 0.015, 0.536,1), [Vec3(0,0,1)], [90]))
    

    #companion cube
    cube_mat = Material(Lambertian(), 0.0, Texture("data/companionCube.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 0.25, Vec4(0, -0.75, 2,1), [Vec3(0,1,0)], [0]))



    # right portal background
    portal_mat = Material(Lambertian(), 0.0,  nothing, RGB{Float32}(0.0, 0.47058823529, 1))
    append!(objs, mesh_helper(portal_background_mesh(2,1), portal_mat, 0.49, Vec4(1.49, 0, 2.0,1), [Vec3(0,1,0),Vec3(1,0,0)], [90,90]))

    # right portal
    portal_mat = Material(Lambertian(), 1.0,  nothing, RGB{Float32}(0.0, 0.0, 0.0))
    append!(objs, mesh_helper(portal_mesh(2,1), portal_mat, 0.47, Vec4(1.46, 0.015, 2.0,1), [Vec3(0,1,0),Vec3(1,0,0)], [90,90]))
    
    lights = [ 
               DirectionalLight(0.8, Vec3(-1,1,1)),
               PointLight(0.9, Vec3(0, 4.5, 2)) ]


   return (Scene(bg, objs, lights), camera_portal2(img_height, img_width))


end



function portalScene3(img_height, img_width)

    bg = black

    objs = []
   
    #make the ground
    push!(objs, Sphere(Vec3(0, -5001, 0), 5000, Material(Lambertian(), 0.03, nothing, RGB{Float32}(0.8, 0.8, 1.0))))
    
    
    #left cube
    cube_mat = Material(Lambertian(), 0.0, Texture("data/wall.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 1.0, Vec4(0, 0, 0, 1)))
    
        
    #right cube
    cube_mat = Material(Lambertian(), 0.0, Texture("data/wall.png", false), white)
    append!(objs, mesh_helper(cube_mesh(), cube_mat, 1.0, Vec4(0, 0, 6, 1)))

    
    #add a bunny
    bunny_mat = Material(Lambertian(), 0.0, nothing, RGB{Float32}(0.8, 0.8, 0.8))
    bunny = read_obj("data/bunny.obj")
    append!(objs, mesh_helper(bunny, bunny_mat, 0.5, Vec4(0, -0.5, 3, 1)))




    #left portal background
    portal_mat = Material(Lambertian(), 0.0,  nothing, RGB{Float32}(0.0, 0.47058823529, 1))
    append!(objs, mesh_helper(portal_background_mesh(2,1), portal_mat, 0.49, Vec4(0, 0, 0.52,1), [Vec3(0,0,1)], [90]))

    #left portal
    portal_mat = Material(Lambertian(), 1.0,  nothing, RGB{Float32}(0.0, 0.0, 0.0))
    append!(objs, mesh_helper(portal_mesh(2,1), portal_mat, 0.47, Vec4(0, 0.015, 0.558,1), [Vec3(0,0,1)], [90]))



    # right portal background
    portal_mat = Material(Lambertian(), 0.0,  nothing, RGB{Float32}(0.99215686274, 0.4, 0.0))
    append!(objs, mesh_helper(portal_background_mesh(2,1), portal_mat, 0.49, Vec4(0, 0, 5.48,1), [Vec3(0,0,1), Vec3(0,1,0)], [90, 180]))


    #right portal
    portal_mat = Material(Lambertian(), 1.0,  nothing, RGB{Float32}(0.0, 0.0, 0.0))
    append!(objs, mesh_helper(portal_mesh(2,1), portal_mat, 0.47, Vec4(0, 0.015, 5.455,1), [Vec3(0,0,1), Vec3(0,1,0)], [90, 180]))

    
    lights = [ 
        PointLight(1.0, Vec3(0,3,3)),
	DirectionalLight(1.0, Vec3(-1, 0.5, 0)) ]


   return (Scene(bg, objs, lights), camera_portal3(img_height, img_width))


end



scenes = [portalScene1, portalScene2, portalScene3]

end # module TestScenes
