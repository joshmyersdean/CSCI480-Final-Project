module OBJMeshes

export read_obj, write_obj, tri_vertex_str
export cube_mesh, cylinder_mesh, sphere_mesh, estimate_normals
export OBJTriangle, OBJMesh

using FileIO
using LinearAlgebra

#push!(LOAD_PATH, pwd())
using ..GfxBase


""" OBJTriangle
A struct that represents a single triangle in a mesh. """
mutable struct OBJTriangle
    positions::Array{Int, 1} # vertex position indices
    uvs::Array{Int, 1} # vertex texture coordinate indices
    normals::Array{Int, 1} # normal vector indices
end

""" OBJMesh
A struct that represents an indexed triangle mesh for reading from
or writing to OBJ format. """
mutable struct OBJMesh
    positions::Array{Vec3, 1} # all vertex positions
    uvs::Array{Vec2, 1} # all texture coordinates
    normals::Array{Vec3, 1} # all vertex normals
    triangles::Array{OBJTriangle, 1} # the OBJTriangles belonging to the mesh
end


""" read_obj(obj_filename)
Read a mesh in OBJ format from file obj_filename."""
function read_obj(obj_filename)
    m = OBJMesh([], [], [], []) # create a mesh
    open(obj_filename) do f
        for (line_number, line) in enumerate(eachline(f))
            if line == "" || line[1] == "#"
                continue # skip comments
            end
            # Read the line and add its contents to the correct field of m:
            tokens = split(strip(line))
            if tokens[1] == "v" # vertex
                push!(m.positions, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vt" # vertex texture
                push!(m.uvs, Vec2([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vn" # vertex normal
                push!(m.normals, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "f"
                # create a OBJTriangle face:
                points = []
                uvs = []
                normals = []
                # handle faces with no texture and/or normals
                for corner in tokens[2:end]
                    indices = split(corner, '/')
                    if length(indices) == 3 # all 3 present, third is a normal
                        push!(normals, parse(Int, indices[3]))
                    end
                    if length(indices) >= 2 && indices[2] != ""
                        # if there are 2 or more and the second isn't blank, it's a texture
                        push!(uvs, parse(Int, indices[2]))
                    end
                    if length(indices) >= 1 # first value is the position
                        push!(points, parse(Int, indices[1]))
                    else # unless it has none, in which case it's not valid
                        error("in line $line_number: face vertex $corner could not be parsed")
                    end
                end
                # create the triangle and add it to the triangles array
                push!(m.triangles, OBJTriangle(points, uvs, normals))
            end
        end
    end
    return m
end

""" write_obj(obj_filename)
Write the given mesh in OBJ format to file obj_filename."""
function write_obj(obj_filename, mesh::OBJMesh)
    open(obj_filename, "w") do f
        # write all positions:
        for v in mesh.positions
            write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all texture coords:
        for v in mesh.uvs
            write(f, "vt $(v[1]) $(v[2])\n")
        end
        # write all normals:
        for v in mesh.normals
            write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
        end

        # write all triangles:
        for tri in mesh.triangles
            write(f, "f $(tri_vertex_str(tri))\n")
        end

    end

end

""" tri_vertex_str(triangle)
Return a string with the indices of applicable positions, texture coordinates,
and normals for a given triangle according to the OBJ specification.
In particular, if p, u, and n are position, vertex and normal, each corner
of the triangle is represented as one of the following:
    p       (position only)
    p/u     (position and texture)
    p//n    (position and normal)
    p/u/n   (position, texture, and normal)
"""
function tri_vertex_str(triangle::OBJTriangle)
    # determine whether textures and normals are present:
    write_uv = length(triangle.uvs) == length(triangle.positions)
    write_normals = length(triangle.normals) == length(triangle.positions)
    corners = []
    for i = 1:3
        output = "$(triangle.positions[i])"
        if write_uv && !write_normals
            output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
        elseif !write_uv && write_normals
            output = output * "//$(triangle.normals[i])"
        elseif write_uv && write_normals
            output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
        end
        push!(corners, output)
    end
    join(corners, " ")
end


""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. """
function cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3( 1, -1, -1)) # 1 RBC
    push!(positions, Vec3( 1, -1,  1)) # 2 RBF
    push!(positions, Vec3(-1, -1,  1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3( 1,  1, -1)) # 5 RTC
    push!(positions, Vec3( 1,  1,  1)) # 6 RTF
    push!(positions, Vec3(-1,  1,  1)) # 7 LTF
    push!(positions, Vec3(-1,  1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3( 1, 0, 0)) # R
    push!(normals, Vec3(-1, 0, 0)) # L
    push!(normals, Vec3( 0, 1, 0)) # U
    push!(normals, Vec3( 0,-1, 0)) # D
    push!(normals, Vec3( 0, 0, 1)) # C
    push!(normals, Vec3( 0, 0,-1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1,2,3], [1,2,3], [4,4,4])) # bottom face 1
    push!(triangles, OBJTriangle([1,3,4], [1,3,4], [4,4,4])) # bottom face 2
    push!(triangles, OBJTriangle([1,5,6], [4,1,2], [1,1,1])) # right face 1
    push!(triangles, OBJTriangle([1,6,2], [4,2,3], [1,1,1])) # right face 2
    push!(triangles, OBJTriangle([2,6,7], [4,1,2], [5,5,5])) # far face 1
    push!(triangles, OBJTriangle([2,7,3], [4,2,3], [5,5,5])) # far face 2
    push!(triangles, OBJTriangle([3,7,8], [2,3,4], [2,2,2])) # left face 1
    push!(triangles, OBJTriangle([3,8,4], [2,4,1], [2,2,2])) # left face 2
    push!(triangles, OBJTriangle([4,8,5], [2,3,4], [6,6,6])) # far face 1
    push!(triangles, OBJTriangle([4,5,1], [2,4,1], [6,6,6])) # far face 2
    push!(triangles, OBJTriangle([5,8,7], [1,2,3], [3,3,3])) # top face 1
    push!(triangles, OBJTriangle([5,7,6], [1,3,4], [3,3,3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end


""" cylinder_mesh(n)
Return a new OBJMesh object approximation of a cylinder with radius 1 and
height 2, centered at the origin. The logitudinal axis is aligned with y, and
it is tesselated with n divisions arranged radially around the outer surface.
The ends of the cylinder are disc-shaped caps parallel to the xz plane.
"""

function portal_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = [] 
    push!(positions, Vec3(0,0,1)) # top of cap
    push!(normals, Vec3(0,0,1))  
    n = 16 
    for i = 1:n+1
        ratio = i/n  
        theta = 2*pi*ratio
        push!(positions, Vec3(-sin(theta)*.75, -cos(theta)*1.25, 1))
        push!(normals, Vec3(-sin(theta),  -cos(theta), 1))
    end
	length_pos = length(positions)
     for i=2:length_pos
       push!(triangles, OBJTriangle([i, 1,i+1],[i+length_pos,1,i+length_pos],[1,1,1]))
     end
     OBJMesh(positions,uvs ,normals,triangles) 
end




function cylinder_mesh(divisionsU)
  positions = []
  uvs = []
  normals = []
  triangles = []

  n = divisionsU # change to n for ease of typing
  push!(positions, Vec3(0,1,0)) # top of cap
  push!(normals, Vec3(0,1,0)) 

  push!(positions, Vec3(0,-1,0)) # bottom of cap
  push!(normals, Vec3(0,-1,0)) 
  
  push!(positions, Vec3(0,1,-1)) # start of shell
  push!(positions, Vec3(0,-1,-1))
  push!(normals, Vec3(0,0,-1)) 

  push!(uvs, Vec2(.75, .75)) # initial uv coordinates
  push!(uvs, Vec2(.25, .75))
  push!(uvs, Vec2(0, .5))
  push!(uvs, Vec2(0, 0))
  # positions
  for i = 1:n-1
    ratio = i/n  
    theta = 2*pi*ratio
    push!(positions, Vec3(-sin(theta), 1, -cos(theta)))
    push!(positions, Vec3(-sin(theta), -1, -cos(theta)))
    push!(normals, Vec3(-sin(theta), 0, -cos(theta)))
    push!(uvs, Vec2(ratio, .5)) # uvs for shell
    push!(uvs, Vec2(ratio, 0))
  end

  length_pos = length(positions) # for later use
  length_norm = length(normals)

  push!(uvs, Vec2(1, .5)) # missing coordinates
  push!(uvs, Vec2(1, 0))
  push!(uvs, Vec2(.75, 1)) 
  push!(uvs, Vec2(.25, .5))

  # uvs for caps
  for i=5:length_pos
      x = positions[i][1] # -sin(θ)
      y = positions[i][2] # y ∈ (-1,1)
      z = positions[i][3] # -cos(θ)
    if y == 1.0 # top cap
      push!(uvs, Vec2((x+1)/4+.5, (1-z)/4+.5))
    else # bottom cap
      push!(uvs, Vec2((x+1)/4, (z+1)/4+.5))
    end
  end

  n1 = 3 # normal counters
  n2 = 3

  # create the shell ignoring the top and bottom coordinates
  for i =3:2:length_pos-3
    # create triangles in a ccw order, texture coords are same index
    push!(triangles, OBJTriangle([i,i+1, i+2],[i,i+1,i+2],
                                 [n1,n1,n1+1]))
    push!(triangles, OBJTriangle([i+1,i+3, i+2],[i+1,i+3,i+2]
                                 ,[n2, n2+1, n2+1])) 
    n1 += 1 # normals increment half as fast as the shell
    n2 += 1
  end
  
  # last bit of the shell
  push!(triangles, OBJTriangle([length_pos-1,length_pos,3],[length_pos-1,length_pos,length_pos+1],
                               [length_norm, length_norm, 3]))
  push!(triangles, OBJTriangle([length_pos,4,3],[length_pos, length_pos+2, length_pos+1],
                               [length_norm,3 ,3 ]))
  
  # create the caps for the cylinder
  for i=3:2:length_pos-2
    push!(triangles, OBJTriangle([i, i+2, 1],[i+length_pos,i+length_pos+2,1],[1,1,1]))
    push!(triangles, OBJTriangle([i+1, 2, i+3],[i+length_pos+1,2,i+length_pos+3],[2,2,2]))
  end

  # last bit of the caps
  push!(triangles, OBJTriangle([length_pos-1, 3, 1], [2*length_pos-1,length_pos+3,1], [1,1,1]))
  push!(triangles, OBJTriangle([length_pos, 2, 4], [2*length_pos,2,length_pos+4], [2,2,2]))

  OBJMesh(positions,uvs ,normals,triangles)
end


""" sphere_mesh(n, m
Create a Latitude-Longitude-tesselated approximation of a sphere with radius 1
centered at the origin. There are n divisions around the equator and m
divisions from pole to pole along each line of longitude. The North pole is at
(0,1,0), the South pole at (0,-1,0)((-1, and points on the Greenwich meridian are
in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
coordinate varies with latitude with v=0 at the South pole and v=1 at the North
pole. Normals should be normal to the ideal sphere surface. See the assignment
for a diagram and further details. """
function sphere_mesh(n, m)
    positions = []
    uvs = []
    normals = []
    triangles = []

    u = 2 * pi / n # normalize 
    v = pi / m # normalize
    
    push!(positions, Vec3(0, -1, 0)) # bottom pole
    push!(normals, Vec3(0,-1,0)) 

    push!(positions, Vec3(0, 1, 0)) # top pole
    push!(normals, Vec3(0,1,0))    
    
    # In this next block we convert x,y,z to spherical coordinates and 
    # add the vectors to the positions array
    for i=1:m-1  
      theta = (pi / 2) - (i * v) # π/2 to -π/2
      y = sin(theta) # Spherical y coord 
      for j=0:n-1
        phi = j * u
        x = cos(theta) * sin(phi) # Spherical x coord
        z = cos(theta) * cos(phi) # spherical y coord
        push!(positions, Vec3(-x, -y, -z)) # negative as we start from the bottom
        push!(normals, Vec3(-x,-y,-z))
      end
    end
    length_pos = length(positions)
    # uv coords
    for i=0:m
      for j=0:n
        push!(uvs, Vec2(j/n,i/m))
      end
    end
   
    # create the bottom "ring"
    for i=1:n-1
      push!(triangles, OBJTriangle([1, i+3,i+2],[i,i+n+1,i+n],[1,i+3,i+2]))
    end

    # missing piece of bottom ring
    push!(triangles, OBJTriangle([1,3,n+2],[n, 2*n+1,2*n],[1,3,n+2]))
    
    t1, t2 = n+1,2*n+2 # texture counters
    for i=3:length_pos-n
      # if i is a multiple of m we subtract m from the first and second coordinates, respectively.
      if ((i-2)%n) == 0
       push!(triangles,OBJTriangle([i,i-n+1,i+n],[t1,t1+1,t2],[i,i-n+1,i+n]))
       push!(triangles,OBJTriangle([i-n+1,i+1,i+n],[t1+1,t2+1,t2],[i-n+1,i+1,i+n]))
       t1 +=2 # textures increment twice as fast
       t2 += 2
       continue
      end

     # for most of the triangles, i is not a multiple of m, add m to 3rd coord.
     push!(triangles, OBJTriangle([i, i+1,i+n],[t1,t1+1,t2],[i, i+1,i+n]))
     push!(triangles, OBJTriangle([i+1,i+1+n,i+n],[t1+1,t2+1,t2],[i+1,i+1+n,i+n]))
     t1 += 1 # textures increment normally here
     t2 += 1
    end
    
    # top "ring"
    for i =length_pos-n+1:length_pos-1
      push!(triangles,OBJTriangle([i,i+1,2],[t1,t1+1,t1+1+n],[i,i+1,2]))
      t1 += 1
      t2 += 1
    end
    # missing bit of the top "ring"
    push!(triangles,OBJTriangle([length_pos,length_pos-n+1,2],[t1,t1+1,t1+1+n],
                                [length_pos, length_pos-n+1,2]))

   deleteat!(uvs, 1) # delete the first/last index of texture coords as it is unneeded
   deleteat!(uvs, length(uvs))
   OBJMesh(positions,uvs, normals, triangles)
end

""" 
    estimate_normals(mesh::OBJMesh)
Estimates normals for the given mesh. Overwrites any existing normals and returns a new OBJMesh object.
"""
function estimate_normals(mesh::OBJMesh)
  normals = []
  # cerate a normal for every position
  for  x in mesh.positions
    push!(normals, Vec3(0,0,0))
  end
  # average all the normals for every triangle
  for x in mesh.triangles
    a = mesh.positions[x.positions[1]]
    b = mesh.positions[x.positions[2]]
    c = mesh.positions[x.positions[3]]
    # definition of computing a normal
    cr = cross(b-a, c-a)
    normals[x.positions[1]] += cr/norm(cr)
    normals[x.positions[2]] += cr/norm(cr)
    normals[x.positions[3]] += cr/norm(cr)
  end

  tmp = [] # create another array to normalize the normals
  for i in normals
    i = i/norm(i) # normalize
    push!(tmp, i)
  end

  normals = tmp # for the return value
  for x in mesh.triangles
    x.normals = x.positions # normals have same index as positions per the reference OBJs
  end
  OBJMesh(mesh.positions, mesh.uvs, normals, mesh.triangles)
end

end # module OBJMeshes
