#!/usr/bin/env julia
module blockster

using ArgParse
import JSON
import DataStructures

using Blockster

function main(args)
    s = ArgParseSettings(description =
            "Blockster. A package for generating multi-block structured hexadral meshes")

    @add_arg_table s begin
         "--dictionary", "-d"        
             help = "Dictionary defining the mesh."
             required = false
             default = joinpath("tests", "cube.json")
         "--nowrite"
             action = :store_true
             help = "Do not write the mesh. For performance tests"
         "--intsize"
             help = "Number of bits used for integers"
             required = false
             arg_type = DataType
             default = Int32
    end

    parsedArgs::Dict{String, Any} = parse_args(s) 


    #const Label::DataType = eval(parse("Int$(parsedArgs["intsize"])"))
    const Label::DataType = parsedArgs["intsize"]

    dictPath::String = parsedArgs["dictionary"]
    nowrite::Bool = parsedArgs["nowrite"]

    # Parse the dictionary
    dict::DataStructures.OrderedDict = JSON.parsefile(dictPath,
                                                      dicttype=DataStructures.OrderedDict)
    # Parse user difined variables
    if haskey(dict, "variables")
        variables = dict["variables"]
        variablesAsStrings = Blockster.variables_as_strings(variables)
    else
        variablesAsStrings = Vector{String}(0)
    end

    # Get the number of blocks
    nBlocks::Label = length(dict["blocks"])
    nVertices::Label = length(dict["vertices"])

    # Vertices defining the mesh as defined in the dict
    vertices = dict["vertices"]
    vertices = Blockster.parse_vertices(variablesAsStrings, vertices) 

    patchNames, patchTypes, patchSurfaces = Blockster.read_boundary(dict, Label)
    patchDicts = dict["boundary"]

    Blockster.check_patch_vertex_labels(patchNames, patchSurfaces, vertices)

    print("Creating blocks...")
    blocks = Blockster.create_blocks(dict["blocks"], vertices, variablesAsStrings, Label)
    print(" Done\n")

    # Create mesh from the blocks
    # blocks as cells
    blocksAsCells = [convert(Blockster.Cell{Label}, blocks[i]) for i in 1:nBlocks]

    # blocks as surface list
    blocksAsSurfaceList = [Blockster.cellfaces(blocksAsCells[i]) for i in 1:nBlocks]

    # Create vertex to block adressing
    vertexBlockAddressing = Blockster.point_cell_addressing(blocksAsCells, nVertices)

    print("Creating block topology...")

    defaultSurfaceStart,
    surfaces,
    blocksAsSurfaces = Blockster.create_topology(
                        blocksAsCells,
                        patchSurfaces,
                        patchNames,
                        vertexBlockAddressing,
                        nVertices
                    )[3:end]
    print(" Done\n")

    nSurfaces::Label = length(surfaces)

    nDefaultSurfaces = nSurfaces - defaultSurfaceStart

    if nDefaultSurfaces > 0
        error("Some boundary surfaces are not part of any patch.")
    end

    ownerBlocks, neighbourBlocks =
        Blockster.create_owner_neighbour(nSurfaces, blocksAsSurfaces)

    nInternalSurfaces::Label = length(neighbourBlocks)

    print("Creating merge list...")
    nCells, nPoints, blockOffsets, mergeList =
        Blockster.calc_merge_info(
            blocks,
            surfaces,
            blocksAsSurfaces,
            ownerBlocks,
            neighbourBlocks,
            nInternalSurfaces
        )
    print(" Done\n")

    print("Creating global point list...")
    points = Blockster.create_points(
                 blocks,
                 blockOffsets,
                 mergeList,
                 nPoints
             )
    print(" Done\n")

    print("Creating global cell list...")
    cells = Blockster.create_cells(
                blocks,
                blockOffsets,
                mergeList,
                nCells
            )
    print(" Done\n")

    print("Creating patches...")
    patches = Blockster.create_patches(
                  blocks,
                  blockOffsets,
                  mergeList,
                  blocksAsSurfaceList,
                  patchSurfaces,
                  surfaces,
                  ownerBlocks
              )
    print(" Done\n")

    pointCellAddressing = Blockster.point_cell_addressing(cells, nPoints)

    print("Creating mesh topology...")

    patchSizes,
    patchStarts,
    defaultPatchStart,
    faces,
    cellsAsFaces = Blockster.create_topology(
                       cells,
                       patches,
                       patchNames,
                       pointCellAddressing,
                       nPoints
                   )
    print(" Done\n")

    nFaces::Label  = length(faces)

    print("Creating owner and neighbour lists...")
    owner, neighbour = Blockster.create_owner_neighbour(nFaces, cellsAsFaces)
    print(" Done\n")

    nInternalFaces = length(neighbour)
    
    # change to 0-based arrays
    owner -= 1
    neighbour -= 1 
    patchStarts -= 1
    faces -= 1

    println("Writing mesh")

    if !isdir(joinpath(".", "test_case", "constant", "polyMesh"))
       mkpath(joinpath(".", "test_case", "constant", "polyMesh"))
    end

    writeDir = joinpath(".", "test_case", "constant", "polyMesh")

    if !nowrite
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
    end

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
