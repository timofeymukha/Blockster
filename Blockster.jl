__precompile__() 
module Blockster

import JSON
import DataStructures
using ArgParse
using StaticArrays: SVector

include("MeshPrimitives.jl")
include("Edges.jl")
include("Blocks.jl")
include("MeshCreate.jl")

#using .MeshCreate
#using .MeshPrimitives
#using .Edges
#using .Blocks


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

function read_boundary!(meshDict)

    nPatches = size(meshDict["boundary"], 1)

    patchNames = Vector{String}(nPatches)
    patchSurfaces = Vector{Vector{Face}}(nPatches)

    for i in 1:nPatches
        patchNames[i] = meshDict["boundary"][i]["name"]
        patchSurfaces[i] = meshDict["boundary"][i]["faces"] + 1
    end

    return patchNames, patchSurfaces

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




function write_header(file, class, location, object)
    write(file, "FoamFile\n")
    write(file, "{\n")
    write(file, "    version     2.0;\n")
    write(file, "    format      ascii;\n")
    write(file, "    class       $class;\n")
    write(file, "    location    \"$location\";\n")
    write(file, "    object      $object;\n")
    write(file, "}\n\n")
end

function create_blocks(dict, vertices)
    nBlocks = size(dict["blocks"], 1)
    blocks = Vector{Block}(nBlocks)

    for blockI in 1:nBlocks
        #println("    Block number ", blockI)
        blocks[blockI] = Block()
        blocks[blockI].nCells = dict["blocks"][blockI]["number of cells"]

        # read vertex numbers defining the block
        blocks[blockI].vertexLabels = dict["blocks"][blockI]["vertices"] + 1

        # the coordinates of all vertices defining the block
        # corresponding to the vertex numbers defining the block
        for j in 1:8
            blocks[blockI].vertices[j] = vertices[blocks[blockI].vertexLabels[j]]
        end

        #println("        Creating edge-points")
        make_block_edges!(blocks[blockI])

        #println("        Creating points")
        create_points!(blocks[blockI])

        #println("        Creating cells")
        create_cells!(blocks[blockI])

        #println("        Creating boundary faces")
        create_boundary_faces!(blocks[blockI])
        #println("        Done")
    end

    return blocks
end

function write_mesh_information(
    nPoints,
    nCells,
    nFaces,
    nInternalFaces,
    patchNames,
    patchStarts,
    patchSizes
)
    println("----------------")
    println("Mesh Information")
    println("----------------")
    println("  nPoints $(nPoints)")
    println("  nCells $(nCells)")
    println("  nFaces $(nFaces)")
    println("  nInternalFaces $(nInternalFaces)")

    println("----------------")
    println("Patches")
    println("----------------")

    for i in 1:length(patchSizes)
        println("  patch $i (start: $(patchStarts[i]) size: $(patchSizes[i])) name: $(patchNames[i])")
    end
end


function parse_vertices(varsAsStr::Vector{String}, vertices)
    n = size(vertices, 1)

    for i in 1:length(varsAsStr)
        eval(parse(varsAsStr[i]))
    end

    floatVertices = Vector{Point}(n)
    for v in 1:n
        currVertex = Vector{Float64}([0, 0, 0])
        for i in 1:3
            if typeof(vertices[v][i]) == String
                currVertex[i] = eval(parse(vertices[v][i]))
            else
                currVertex[i] = vertices[v][i]
            end
        end
        floatVertices[v] = currVertex
    end

    return floatVertices
end

function variables_as_strings(variables)

    varsAsStr = Vector{String}(length(variables))

    i = 1
    for (key, value) in variables
        varsAsStr[i] = "$key = $value"
        i += 1
    end

    return varsAsStr
end

function main(args)
    s = ArgParseSettings(description =
            "Blockster. A package for generating multi-block structured hexadral meshes")

    @add_arg_table s begin
         "--dictionary", "-d"        
             help = "Dictionary defining the mesh."
             required = false
             default = joinpath("tests", "cube.json")
    end

    parsedArgs = parse_args(s) # the result is a Dict{String,Any}

    dictPath = parsedArgs["dictionary"]

    # Parse the dictionary
    dict = JSON.parsefile(dictPath, dicttype=DataStructures.OrderedDict)
    
    # Parse user difined variables
    variables = dict["variables"]
    variablesAsStrings = variables_as_strings(variables)
    println(variablesAsStrings)

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

    patchNames, patchSurfaces =  read_boundary!(dict)
    check_patch_vertex_labels(patchNames, patchSurfaces, vertices)

    println("Creating blocks")
    blocks = create_blocks(dict, vertices)

    # Create mesh from the blocks
    # blocks as cells
    blocksAsCells = [convert(Cell, blocks[i]) for i in 1:nBlocks]

    # blocks as faces 
    blockFaces = [cellfaces(blocksAsCells[i]) for i in 1:nBlocks]

    # Create vertex to block adressing
    vertexBlockAddressing = point_cell_addressing(blocksAsCells, nVertices)



    patchSizes, patchStarts, defaultPatchStart, faces, nFaces, cellsAsFaces =
        create_topology(blocksAsCells,
                        patchSurfaces,
                        patchNames,
                        vertexBlockAddressing,
                        nVertices)

    nDefaultFaces = nFaces - defaultPatchStart

    if nDefaultFaces > 0
        warn("Undefined block faces present in the mesh description")
    end

    owner, neighbour, nInternalFaces = init_mesh(faces, cellsAsFaces)


    nCells, nPoints, blockOffsets, mergeList =
        calc_merge_info(blocks, vertices, faces, cellsAsFaces, owner, neighbour,
                   nInternalFaces)
    points = create_points(blocks, blockOffsets, mergeList, nPoints)
    cells = create_cells(blocks, blockOffsets, mergeList, nCells)
    patches = create_patches(blocks, blockOffsets, mergeList, blockFaces,
                             patchSurfaces, faces, owner)

    pointCellAddressing = point_cell_addressing(cells, nPoints)
    patchSizes, patchStarts, defaultPatchStart, faces, nFaces, cellsAsFaces = 
    create_topology(cells,
                    patches,
                    patchNames,
                    pointCellAddressing,
                    size(points, 1))
    owner, neighbour, nInternalFaces = init_mesh(faces, cellsAsFaces)


    println("Writing mesh")
    # change to 0-based arrays
    owner -= 1
    neighbour -= 1 
    for i in 1:size(faces, 1)
        faces[i] = faces[i] - 1
    end


    if !isdir(joinpath(".", "test_case", "constant", "polyMesh"))
       mkpath(joinpath(".", "test_case", "constant", "polyMesh"))
    end

    ownerFile = open(joinpath(".", "test_case", "constant", "polyMesh", "owner"), "w")
    write_header(ownerFile, "labelList", "constant/polyMesh", "owner")
    write(ownerFile, "$(size(owner, 1))\n") 
    write(ownerFile, "(\n")
    for i in 1:size(owner, 1)
        write(ownerFile, "$(owner[i])\n")
    end
    write(ownerFile, ")\n")
    close(ownerFile)

    neighbourFile = open(joinpath(".", "test_case", "constant", "polyMesh", "neighbour"), "w")
    write_header(neighbourFile, "labelList", "constant/polyMesh", "neighbour")
    write(neighbourFile, "$(size(neighbour, 1))\n") 
    write(neighbourFile, "(\n")
    for i in 1:size(neighbour, 1)
        write(neighbourFile, "$(neighbour[i])\n")
    end
    write(neighbourFile, ")\n")
    close(neighbourFile)

    pointsFile = open(joinpath(".", "test_case", "constant", "polyMesh", "points"), "w")
    write_header(pointsFile, "vectorField", "constant/polyMesh", "points")
    write(pointsFile, "$(size(points, 1))\n") 
    write(pointsFile, "(\n")
    for i in 1:size(points, 1)
        write(pointsFile, "($(points[i][1]) $(points[i][2]) $(points[i][3]))\n")
    end
    write(pointsFile, ")\n")
    close(pointsFile)

    facesFile = open(joinpath(".", "test_case", "constant", "polyMesh", "faces"), "w")
    write_header(facesFile, "faceList", "constant/polyMesh", "faces")
    write(facesFile, "$(size(faces, 1))\n") 
    write(facesFile, "(\n")
    for i in 1:size(faces, 1)
        write(facesFile, "4($(faces[i][1]) $(faces[i][2]) $(faces[i][3]) $(faces[i][4]))\n")
    end
    write(facesFile, ")\n")
    close(facesFile)

    boundaryFile = open(joinpath(".", "test_case", "constant", "polyMesh", "boundary"), "w")
    write_header(boundaryFile, "polyBoundaryMesh", "constant/polyMesh", "boundary")
    write(boundaryFile, "$(size(patches, 1))\n") 
    write(boundaryFile, "(\n")
    for i in 1:size(patches, 1)
        write(boundaryFile, "    $(patchNames[i])\n")
        write(boundaryFile, "    {\n")
        write(boundaryFile, "    type      wall;\n")
#        write(boundaryFile, "    inGroups\n")
        write(boundaryFile, "    nFaces    $(patchSizes[i]);\n")
        write(boundaryFile, "    startFace $(patchStarts[i] - 1);\n")
        write(boundaryFile, "    }\n")
    end

    write(boundaryFile, ")\n")
    close(boundaryFile)

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
