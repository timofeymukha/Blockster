export read_boundary, parse_vertices, parse_ncells, variables_as_strings

function read_boundary(meshDict)

    nPatches = size(meshDict["boundary"], 1)

    patchNames = Vector{String}(nPatches)
    patchSurfaces = Vector{Vector{Face}}(nPatches)

    for i in 1:nPatches
        patchNames[i] = meshDict["boundary"][i]["name"]
        patchSurfaces[i] = meshDict["boundary"][i]["faces"] + 1
    end

    return patchNames, patchSurfaces

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

function parse_ncells(varsAsStr::Vector{String}, nCells)
    n = 3

    for i in 1:length(varsAsStr)
        eval(parse(varsAsStr[i]))
    end

    intNCells = Vector{Int32}(n)
    for dir in 1:n
        if typeof(nCells[dir]) == String
            intNCells[dir]  = eval(parse(nCells[dir]))
        else
            intNCells[dir] = nCells[dir]
        end
    end

    return intNCells
end

function parse_grading(varsAsStr::Vector{String}, grading, gradingType)
    n = 12

    for i in 1:length(varsAsStr)
        eval(parse(varsAsStr[i]))
    end

    # convert grading to full edge grading if necessary
    edgeGrading = Vector{Any}(n)
    
    if gradingType == "simple"
        [edgeGrading[i] =grading[1] for i in 1:4]
        [edgeGrading[i] =grading[2] for i in 5:8]
        [edgeGrading[i] =grading[3] for i in 9:12]
    else
        edgeGrading = grading
    end

    for i in 1:n
        if !(typeof(edgeGrading[i]) <: AbstractArray)
            edgeGrading[i] = Vector([Vector([1., 1., edgeGrading[i]])])
        end
    end

    for edge in 1:n
        for section in 1:length(edgeGrading[edge])
            for i in 1:3
                current = @view edgeGrading[edge][section][i]
                if typeof(current[1]) == String
                    current[1]  = Float64(eval(parse(current[1])))
                else
                    current[1]  = Float64(current[1])
                end
            end
        end
    end

    return edgeGrading
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
