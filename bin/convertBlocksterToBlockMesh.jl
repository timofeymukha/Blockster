module ToBlockMesh

import JSON
import DataStructures
using ArgParse

include("../src/Blockster.jl")

function main(args)
    s = ArgParseSettings(description =
            "Convert a Blockster mesh specification dictionary to a blockMeshDict.")

    @add_arg_table s begin
         "--dictionary", "-d"        
             help = "Dictionary defining the mesh."
             required = false
             default = joinpath("tests", "dict.json")
    end

    parsedArgs = parse_args(s) 

    dictPath = parsedArgs["dictionary"]
    # Parse the dictionary
    dict = JSON.parsefile(dictPath, dicttype=DataStructures.OrderedDict)
    
    # Parse user difined variables
    if haskey(dict, "variables")
        variables = dict["variables"]
        variablesAsStrings = Blockster.variables_as_strings(variables)
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

    blockMeshDict = open(joinpath("./test_case", "system", "blockMeshDict"), "w")
    Blockster.write_header(blockMeshDict, "dictionary", "system", "blockMeshDict")

    write(blockMeshDict, "vertices\n(\n")
    for i in 1:length(vertices)
        write(blockMeshDict, "  ($(vertices[i][1]) $(vertices[i][2]) $(vertices[i][3]))\n")
    end
    write(blockMeshDict, ");\n\n")


    write(blockMeshDict, "edges\n(\n")
    write(blockMeshDict, ");\n\n")

    write(blockMeshDict, "blocks\n(\n")
    for i in 1:length(dict["blocks"])
        nCells = Blockster.parse_ncells(variablesAsStrings, dict["blocks"][i][2])
        gradingType = dict["blocks"][i][3]
        grading = Blockster.parse_grading(
                      variablesAsStrings,
                      dict["blocks"][i][4],
                      gradingType
                  )
        write(blockMeshDict, "  hex (")
        for j in 1:8
            write(blockMeshDict, "$(dict["blocks"][i][1][j]) ")
        end
        write(blockMeshDict, ") ($(nCells[1]) $(nCells[2]) $(nCells[3]))\n")

        if gradingType == "simple"
            write(blockMeshDict, "  simpleGrading\n  (\n")
            for dir = (1, 5, 9)
                if length(grading[dir]) == 1
                    write(blockMeshDict, "    $(grading[dir][1][3])\n")
                else 
                    write(blockMeshDict, "    (\n")
                    for secI = 1:length(grading[dir])
                        grI = grading[dir][secI]
                        write(blockMeshDict, "      ($(grI[1]) $(grI[2]) $(grI[3]))\n")

                    end
                    write(blockMeshDict, "    )\n")
                end

            end
            write(blockMeshDict, "  )\n")
        end


    end

    write(blockMeshDict, ");\n\n")

    write_boundary(blockMeshDict, patchDicts)
        

    close(blockMeshDict)

end


function write_boundary(fileName, patchDicts)
    write(fileName, "boundary\n(\n")

    for i in 1:length(patchDicts)
        write(fileName, "  $(patchDicts[i]["name"])\n")
        write(fileName, "  {\n")
        write(fileName, "    type $(patchDicts[i]["type"]);\n")
        for (key, val) in patchDicts[i]
            if !(key in ["name", "type", "faces"])
                write(fileName, "    $(key) $(val);\n")
            end
        end
        write(fileName, "    faces (")

        for fI in 1:length(patchDicts[i]["faces"])
            face = patchDicts[i]["faces"][fI]
            write(fileName, "($(face[1]) $(face[2]) $(face[3]) $(face[4])) ")
        end
        write(fileName, ");\n")
        write(fileName, "  }\n")
    end

    write(fileName, ");\n\n")
end

main(ARGS)
end
