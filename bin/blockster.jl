module blockster

using ArgParse
import JSON
import DataStructures

using Blockster
#include("../src/Blockster.jl")

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
             default = 32
    end

    parsedArgs::Dict{String, Any} = parse_args(s) # the result is a Dict{String,Any}


    const Label::DataType = eval(parse("Int$(parsedArgs["intsize"])"))

#    run(parsedArgs, Label)
    

#end

#function run(parsedArgs::Dict{String, Any}, ::Type{Label}) where {Label <: Integer}
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
    blocks = Blockster.create_blocks(dict, vertices, variablesAsStrings, Label)
    print(" Done\n")

    # Create mesh from the blocks
    # blocks as cells
    blocksAsCells = [convert(Blockster.Cell{Label}, blocks[i]) for i in 1:nBlocks]

    # blocks as faces 
    blockFaces = [Blockster.cellfaces(blocksAsCells[i]) for i in 1:nBlocks]

    # Create vertex to block adressing
    vertexBlockAddressing = Blockster.point_cell_addressing(blocksAsCells, nVertices)

    print("Creating block topology...")
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
    print(" Done\n")

    nDefaultFaces = nFaces - defaultPatchStart

    if nDefaultFaces > 0
        warn("Undefined block faces present in the mesh description")
    end

    owner, neighbour, nInternalFaces = Blockster.init_mesh(faces, cellsAsFaces)


    print("Creating merge list...")
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
    print(" Done\n")

    print("Creating global point list...")
    points = Blockster.create_points(
                 blocks,
                 blockOffsets,
                 mergeList,
                 nPoints
             )
    print(" Done\n")
    nPoints::Label = length(points)

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
                  blockFaces,
                  patchSurfaces,
                  faces,
                  owner
              )
    print(" Done\n")

    pointCellAddressing = Blockster.point_cell_addressing(cells, nPoints)

    print("Creating mesh topology...")

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
                       nPoints
                   )
    print(" Done\n")

    print("Creating owner and neighbour lists...")
    owner, neighbour, nInternalFaces = Blockster.init_mesh(faces, cellsAsFaces)
    print(" Done\n")
    
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
