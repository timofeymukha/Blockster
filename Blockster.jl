__precompile__() 
module Blockster

import JSON
import DataStructures
using ArgParse
using StaticArrays: SVector

include("MeshPrimitives.jl")
include("Read.jl")
include("Edges.jl")
include("Blocks.jl")
include("MeshCreate.jl")
include("Write.jl")

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


function main(args)
    s = ArgParseSettings(description =
            "Blockster. A package for generating multi-block structured hexadral meshes")

    @add_arg_table s begin
         "--dictionary", "-d"        
             help = "Dictionary defining the mesh."
             required = false
             default = joinpath("tests", "channel.json")
    end

    parsedArgs = parse_args(s) # the result is a Dict{String,Any}

    dictPath = parsedArgs["dictionary"]

    # Parse the dictionary
    dict = JSON.parsefile(dictPath, dicttype=DataStructures.OrderedDict)
    
    # Parse user difined variables
    if haskey(dict, "variables")
        variables = dict["variables"]
        variablesAsStrings = variables_as_strings(variables)
    else
        variablesAsStrings = Vector{String}(0)
    end

    # Get the number of blocks
    nBlocks = size(dict["blocks"], 1)
    nVertices = size(dict["vertices"], 1)

    # Vertices defining the mesh as defined in the dict
    #vertices = Vector{Point}(nVertices)
    vertices = dict["vertices"]
    vertices = parse_vertices(variablesAsStrings, vertices) 

    patchNames = String[]
    # The faces of the patches defined from vertices in the dict
    patchSurfaces = Vector{Vector{Face}}(0)

    patchNames, patchSurfaces = read_boundary(dict)
    check_patch_vertex_labels(
        patchNames,
        patchSurfaces,
        vertices
    )

    println("Creating blocks")
    blocks = create_blocks(dict, vertices, variablesAsStrings)

    # Create mesh from the blocks
    # blocks as cells
    blocksAsCells = [convert(Cell, blocks[i]) for i in 1:nBlocks]

    # blocks as faces 
    blockFaces = [cellfaces(blocksAsCells[i]) for i in 1:nBlocks]

    # Create vertex to block adressing
    vertexBlockAddressing = point_cell_addressing(blocksAsCells, nVertices)

    patchSizes,
    patchStarts,
    defaultPatchStart,
    faces,
    nFaces,
    cellsAsFaces = create_topology(
                       blocksAsCells,
                       patchSurfaces,
                       patchNames,
                       vertexBlockAddressing,
                       nVertices
                   )

    nDefaultFaces = nFaces - defaultPatchStart

    if nDefaultFaces > 0
        warn("Undefined block faces present in the mesh description")
    end

    owner, neighbour, nInternalFaces = init_mesh(faces, cellsAsFaces)


    nCells, nPoints, blockOffsets, mergeList =
        calc_merge_info(blocks, vertices, faces, cellsAsFaces, owner, neighbour,
                   nInternalFaces)

    points = create_points(
                 blocks,
                 blockOffsets,
                 mergeList,
                 nPoints
             )

    cells = create_cells(
                blocks,
                blockOffsets,
                mergeList,
                nCells
            )

    patches = create_patches(
                  blocks,
                  blockOffsets,
                  mergeList,
                  blockFaces,
                  patchSurfaces,
                  faces,
                  owner
              )

    pointCellAddressing = point_cell_addressing(cells, nPoints)

    patchSizes,
    patchStarts,
    defaultPatchStart,
    faces,
    nFaces,
    cellsAsFaces = create_topology(
                       cells,
                       patches,
                       patchNames,
                       pointCellAddressing,
                       size(points, 1)
                   )

    owner, neighbour, nInternalFaces = init_mesh(faces, cellsAsFaces)
    
    # change to 0-based arrays
    owner -= 1
    neighbour -= 1 
    for i in 1:size(faces, 1)
        faces[i] = faces[i] - 1
    end

    println("Writing mesh")

    if !isdir(joinpath(".", "test_case", "constant", "polyMesh"))
       mkpath(joinpath(".", "test_case", "constant", "polyMesh"))
    end

    writeDir = joinpath(".", "test_case", "constant", "polyMesh")

    write_mesh(
        writeDir,
        owner,
        neighbour,
        points,
        faces,
        patchNames,
        patchStarts,
        patchSizes
    )

    write_mesh_information(
        nPoints,
        nCells,
        nFaces,
        nInternalFaces, 
        patchNames,
        patchStarts,
        patchSizes
    )

end

main(ARGS)
end
