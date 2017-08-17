__precompile__() 
module Blockster

using StaticArrays: SVector

include("MeshPrimitives.jl")
include("Read.jl")
include("Edges.jl")
include("Blocks.jl")
include("MeshCreate.jl")
include("Write.jl")

"""
    bounding_box(points)

    Compute the bounding box of a set of points.
"""
function bounding_box(points::Vector{Point})
    min = Vector(points[1])
    max = Vector(points[1])

    for i in 2:length(points)
        if points[i][1] < min[1]
            min[1] = points[i][1]
        end
        if points[i][2] < min[2]
            min[2] = points[i][2]
        end
        if points[i][3] < min[3]
            min[3] = points[i][3]
        end
        if points[i][1] > max[1]
            max[1] = points[i][1]
        end
        if points[i][2] > max[2]
            max[2] = points[i][2]
        end
        if points[i][3] > max[3]
            max[3] = points[i][3]
        end
    end
   return min, max 
end


"""
    check_patch_vertex_labels(patchNames, patchSurfaces, vertices)

    Check that the boundaries are defined through vertices that exist.
"""
function check_patch_vertex_labels(patchNames,
                                   patchSurfaces,
                                   vertices)

    for patchI in 1:size(patchNames, 1)
        for faceI in size(patchSurfaces[patchI], 1)
            if !isempty(find(Array(patchSurfaces[patchI][faceI]) .< 1))
                error("""check_patch_vertex_labels() : face vertex label < 0
                       in patch """, patchNames[patchI])
            elseif !isempty(find(Array(patchSurfaces[patchI][faceI]) .>
                                 size(vertices, 1)))
                error("""check_patch_vertex_labels() : face vertex label
                out of bounds in patch """, patchNames[patchI])
            end
        end
    end

end

end
