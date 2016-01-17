module blockMesh

using JSON
using MeshPrimitives
using Edges
using Blocks

# Parse the dictionary
dict = JSON.parsefile("dict.json")

# Get the number of blocks
nBlocks = size(dict["blocks"], 1)
nVertices = size(dict["vertices"], 1)

vertices = Vector{Point}(nVertices)

for i in 1:nVertices
    vertices[i] = dict["vertices"][i]
end


function read_boundary!(meshDict, 
                        patchNames::Array{ASCIIString, 1},
                        patchBlockFaces::Vector{Vector{Face}})

    nPatches = size(meshDict["boundary"], 1)

    resize!(patchNames, nPatches)
    resize!(patchBlockFaces, nPatches)

    # Grab the patch names
    for i in 1:nPatches
        patchNames[i] = meshDict["boundary"][i]["name"]
        patchBlockFaces[i] = meshDict["boundary"][i]["faces"] + 1
    end

end

function check_patch_vertex_labels(patchNames,
                                   patchBlockFaces,
                                   vertices)
    for patchI in 1:size(patchNames, 1)
        for faceI in size(patchBlockFaces[patchI], 1)
            if !isempty(find(Array(patchBlockFaces[patchI][faceI]) .< 1))
                error("""check_patch_vertex_labels() : face vertex label < 0
                       in patch """, patchNames[patchI])
            elseif !isempty(find(Array(patchBlockFaces[patchI][faceI]) .>
                                 size(vertices, 1)))
                error("""check_patch_vertex_labels() : face vertex label 
                out of bounds in patch """, patchNames[patchI])
            end
        end
    end
end

patchNames = ASCIIString[]
patchBlockFaces = Vector{Vector{Face}}(0)

read_boundary!(dict, patchNames, patchBlockFaces)
check_patch_vertex_labels(patchNames, patchBlockFaces, vertices)






"Type defining the mesh"
type Mesh
    "The json file to read the mesh properties from"
    dict

    "The blocks defining the mesh" 
    blocks::Vector{Block}

    "The points of the mesh"
    points::Vector{Point}

    "The vertices of the mesh, used to define the blocks"
    vertices::Vector{Point}
end

println("Creating blocks...")

blocks = Vector{Block}(nBlocks)

function point_cell_addressing(cells::Vector{Cell}, nPoints::Int64)

    # Number of cells each point is included in
    pointCellAddressing = Vector{Vector{Int64}}(nPoints)

    for i in 1:nPoints
        pointCellAddressing[i] = Vector{Int64}(0)
    end

    # For each cells
    for cellI in 1:size(cells, 1)

        # For each point in the cell
        for pointI in 1:size(cells[cellI], 1)

            # Add current cell to the addressing list
            push!(pointCellAddressing[pointI], cellI)
        end

    end

    return pointCellAddressing
end

function patch_face_cells(faces::Vector{Face},
                          cellFaces::Vector{Vector{Face}},
                          pointCellAddressing::Vector{Vector{Int64}})

    faceCells = Vector{Int64}(size(faces,1))

    # For each face determine the cell
    for faceI in 1:size(faces, 1)

        found = false
        currFace = faces[faceI]

        # For each point of the face
        for facePointI in 1:size(currFace, 1)
            
            facePointCells = pointCellAddressing[currFace[facePointI]]

            # For each cell that this point is part of
            for cellI in 1:size(facePointCells, 1)

                cellIFaces = cellFaces[cellI]

                # For each face of the cells check if it is the same as faceI
                for cellIFaceI in 1:size(cellIFaces, 1)

                    if samePoints(cellIFaces[cellIFaceI], faces[faceI])
                        found = true
                        faceCells[faceI] = facePointCells[cellI]
                    end

                    # We got the cell, we can exit all the inner loops
                    if found
                        break
                    end
                end
                if found
                    break
                end
            end
            if found
                break
            end
        end

        if !found
            error("patch_face_cells: Cound not find internal cell.")
        end
    end

    return faceCells
end

function create_topology(cells::Vector{Cell},
                         boundaryFaces::Vector{Vector{Face}},
                         boundaryPatchNames::Vector{ASCIIString},
                         nPoints)

    # Define a vector of cells defined as face vectors
    cellsAsFaces = Vector{Vector{Int64}}(size(cells, 1))

    # Get the faces of each cell and the maximum number of faces
    cellFaces = Vector{Vector{Face}}(size(cells, 1))
    maxFaces = 0

    for i in 1:size(cellFaces,1)

        cellFaces[i] = MeshPrimitives.faces(cells[i])
        maxFaces += size(cellFaces[i],1)

        # Set the associated face label to -1, to mark as undefined
        cellsAsFaces[i] = fill(-1, size(cellFaces[i], 1))
    end

    # Declare the array of faces
    faces = Vector{Face}(maxFaces)

    nFaces = 0

    # Get point to cell addressing
    pointCellAddr = point_cell_addressing(cells, nPoints)

    found = false

    # Go through all cells
    for cellI in 1:size(cells, 1)

        cellIFaces = cellFaces[cellI]

        neiCells = fill(-1, size(cellIFaces))
        faceOfNeiCells = fill(-1, size(cellIFaces))

        nNeighbours = 0
        
        # Identify the neighbours for each face of the cell
        for faceI in 1:size(cellIFaces, 1)
           
            # Skip faces that have already been matched
            if cellsAsFaces[cellI][faceI] > -1
                continue
            end

            found = false

            currFace = cellIFaces[faceI]

            for pointI in 1:size(currFace, 1)

                # Get list of cells sharing this point
                currNeighbours = pointCellAddr[currFace[pointI]]

                for neiI in 1:size(currNeighbours, 1)

                    currNei = currNeighbours[neiI]

                    # Reject neighbours with lower index
                    if currNei < cellI

                        # The list of faces to search through
                        neiFaces = cellFaces[currNei]

                        for neiFaceI in 1:size(neiFaces, 1)
                            
                            if neiFaces[neiFaceI] == currFace
                                #Match!
                                found = true

                                neiCells[faceI] = currNei
                                faceOfNeiCell[faceI] = neiFaceI
                                nNeighbours += 1

                                break
                            end
                        end

                        if found
                            break
                        end
                    end

                    if found
                        break
                    end
                end
                if found
                    break
                end

            end
        end

        # Add the faces in the increasing order of neighbours
        for neiSearch in 1:nNeighbours
            
            # Find the lowest neighbour which is still valid
            nextNei = -1
            minNei = size(cells, 1)

            for ncI in 1:size(neiCells)
                if neiCells[ncI] > -1 && neiCells[ncI] < minNei
                    nextNei = ncI
                    minNei = neiCells[ncI]
                end
            end

            if nextNei > -1
                # Note that nFaces acts as the current face to fill also
                faces[nFaces] = curFaces[nextNei]

                # Set cell-face and cell-neighbour-face to current face label
                cellsAsFaces[cellI][nextNei] = nFaces
                cellsAsFaces[neiCells[nextNei]][faceOfNeiCell[nextNei]] = 
                                                                        nFaces

                # Stop the neighbour from being used again
                neiCells[nextNei] = -1

                # Increment number of faces counter
                nFaces += 1
            else
                error("set_topology(): Error in internal face insertion.")
            end
        end

    end

end

for blockI in 1:nBlocks
    println("    Block number ", blockI)
    blocks[blockI] = Block()
    blocks[blockI].nCells = dict["blocks"][blockI]["number of cells"]
    blocks[blockI].vertexLabels = dict["blocks"][blockI]["vertices"] + 1

    for j in 1:8
        blocks[blockI].vertices[j] = vertices[blocks[blockI].vertexLabels[j]]
    end
    
    # Create the edge points
    println("        Creating edge-points")
    make_block_edges!(blocks[blockI])

    # Create the points
    println("        Creating points")
    create_points!(blocks[blockI])

    # Create the cells
    println("        Creating cells")
    create_cells!(blocks[blockI])
    println("        Done")
end

pointCellAddressing = point_cell_addressing(convert(Vector{Cell}, blocks), size(vertices, 1))
blockFaces = [faces(convert(Cell, blocks[i])) for i in 1:nBlocks]
patchFaceCells = patch_face_cells(patchBlockFaces[1], blockFaces, pointCellAddressing)

blocksAsCells =Vector{Cell}(nBlocks)

for i in 1:nBlocks
    blocksAsCells[i] = convert(Cell, blocks[i])
end

create_topology(blocksAsCells, patchBlockFaces, patchNames, size(vertices, 1)) 


end
