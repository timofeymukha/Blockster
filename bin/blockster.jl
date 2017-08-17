module blockster

using ArgParse
import JSON
import DataStructures

include("../src/Blockster.jl")


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
    vertices = dict["vertices"]
    vertices = Blockster.parse_vertices(variablesAsStrings, vertices) 

    patchNames, patchTypes, patchSurfaces = Blockster.read_boundary(dict)
    patchDicts = dict["boundary"]

    Blockster.check_patch_vertex_labels(patchNames, patchSurfaces, vertices)

    println("Creating blocks")
    blocks = Blockster.create_blocks(dict, vertices, variablesAsStrings)

    # Create mesh from the blocks
    # blocks as cells
    blocksAsCells = [convert(Blockster.Cell, blocks[i]) for i in 1:nBlocks]

    # blocks as faces 
    blockFaces = [Blockster.cellfaces(blocksAsCells[i]) for i in 1:nBlocks]

    # Create vertex to block adressing
    vertexBlockAddressing = Blockster.point_cell_addressing(blocksAsCells, nVertices)

    patchSizes,
    patchStarts,
    defaultPatchStart,
    faces,
    nFaces,
    cellsAsFaces = Blockster.create_topology(
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

    owner, neighbour, nInternalFaces = Blockster.init_mesh(faces, cellsAsFaces)


    nCells, nPoints, blockOffsets, mergeList =
        Blockster.calc_merge_info(
            blocks,
            vertices,
            faces,
            cellsAsFaces,
            owner,
            neighbour,
            nInternalFaces
        )

    points = Blockster.create_points(
                 blocks,
                 blockOffsets,
                 mergeList,
                 nPoints
             )

    cells = Blockster.create_cells(
                blocks,
                blockOffsets,
                mergeList,
                nCells
            )

    patches = Blockster.create_patches(
                  blocks,
                  blockOffsets,
                  mergeList,
                  blockFaces,
                  patchSurfaces,
                  faces,
                  owner
              )

    pointCellAddressing = Blockster.point_cell_addressing(cells, nPoints)

    patchSizes,
    patchStarts,
    defaultPatchStart,
    faces,
    nFaces,
    cellsAsFaces = Blockster.create_topology(
                       cells,
                       patches,
                       patchNames,
                       pointCellAddressing,
                       size(points, 1)
                   )

    owner, neighbour, nInternalFaces = Blockster.init_mesh(faces, cellsAsFaces)
    
    # change to 0-based arrays
    owner -= 1
    neighbour -= 1 
    for i in 1:size(faces, 1)
        faces[i] = faces[i] - 1
    end

    for i in 1:size(patchStarts, 1)
        patchStarts[i] = patchStarts[i] - 1
    end

    println("Writing mesh")

    if !isdir(joinpath(".", "test_case", "constant", "polyMesh"))
       mkpath(joinpath(".", "test_case", "constant", "polyMesh"))
    end

    writeDir = joinpath(".", "test_case", "constant", "polyMesh")

    Blockster.write_mesh(
        writeDir,
        owner,
        neighbour,
        points,
        faces,
        patchStarts,
        patchSizes,
        patchDicts
    )

    Blockster.write_mesh_information(
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
