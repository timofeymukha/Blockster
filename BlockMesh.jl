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
            
            currPoint = cells[cellI][pointI]
            # Add current cell to the addressing list
            push!(pointCellAddressing[currPoint], cellI)
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

            #println(facePointCells)

            # For each cell that this point is part of
            for cellI in 1:size(facePointCells, 1)

                cellIFaces = cellFaces[facePointCells[cellI]]
                #println(cellIFaces)

                # For each face of the cells check if it is the same as faceI
                for cellIFaceI in 1:size(cellIFaces, 1)

                    if samepoints(cellIFaces[cellIFaceI], faces[faceI])
                #       println("match!")
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
            error("patch_face_cells(): Cound not find internal cell for face
                   $currFace.")
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

        cellFaces[i] = cellfaces(cells[i])
        maxFaces += size(cellFaces[i],1)

        # Set the associated face label to -1, to mark as undefined
        cellsAsFaces[i] = fill(-1, size(cellFaces[i], 1))
    end

    # Declare the array of faces
    faces = Vector{Face}(maxFaces)

    nFaces = 1

    # Get point to cell addressing
    pointCellAddr = point_cell_addressing(cells, nPoints)

    #println("nPoints: $nPoints")
    #println("cells: $cells")
    #println("pointCellsAddr: $pointCellAddr")

    found = false

    # Add the non-boundary faces of all cells to the face list
    # Form the cellsAsFaces array
    for cellI in 1:size(cells, 1)

        cellIFaces = cellFaces[cellI]

        neiCells = fill(-1, size(cellIFaces))
        faceOfNeiCell = fill(-1, size(cellIFaces))

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
                    if currNei > cellI

                        # The list of faces to search through
                        neiFaces = cellFaces[currNei]

                        for neiFaceI in 1:size(neiFaces, 1)
                            
                            if samepoints(neiFaces[neiFaceI], currFace)
                            #if neiFaces[neiFaceI] == currFace
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
            minNei = size(cellsAsFaces, 1)

            for ncI in 1:size(neiCells, 1)
                if neiCells[ncI] > -1 && neiCells[ncI] <= minNei
                    nextNei = ncI
                    minNei = neiCells[ncI]
                end
            end

            if nextNei > -1
                # Note that nFaces acts as the current face to fill also
                faces[nFaces] = cellIFaces[nextNei]

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

    # Do the boundary faces

    patchSizes = Vector{Int64}(size(boundaryFaces, 1))
    patchStarts = Vector{Int64}(size(boundaryFaces, 1))

    for patchI in 1:size(boundaryFaces, 1)

        patchFaces = boundaryFaces[patchI]
        patchName = boundaryPatchNames[patchI]

        currPatchFaceCells = patch_face_cells(patchFaces, cellFaces,
                                              pointCellAddressing)

        currPatchStart = nFaces


        for faceI in 1:size(patchFaces, 1)
            
            currFace = patchFaces[faceI]


            # Get the cell on the inside corresponding to this face
            cellInside = currPatchFaceCells[faceI]

            # Get the faces of that cell
            facesOfCellInside = cellFaces[cellInside]

            for cellFaceI = 1:6
                if samepoints(facesOfCellInside[cellFaceI], currFace)
                    if cellsAsFaces[cellInside][cellFaceI] >= 0
                        error("""set_topology(): Trying to specify a boundary face
                              $currFace
                              on the face of cell $cellInside, which is either 
                              an internal face or already belongs to some other
                              patch.
                              
                              This is face $faceI of patch $patchI named
                              $patchName.""")
                    end

                    found = true


                    # Set the patch face to the corresponding cell-face
                    faces[nFaces] = facesOfCellInside[cellFaceI]
                    cellsAsFaces[cellInside][cellFaceI] = nFaces

                    break

                end
            end

            if !found
                error("set_topology(): face does not seem to belong to a cell,
                        which, according to addressing, should be next to it")
            end

            nFaces += 1
        end

        patchSizes[patchI] = nFaces - currPatchStart
        patchStarts[patchI] = currPatchStart

    end

    defaultPatchStart = nFaces

    # Take care of "non-existing faces", put them into the default patch
    for cellI in 1:size(cellsAsFaces, 1)

        currCellFaces = cellsAsFaces[cellI]

        for faceI in 1:size(currCellFaces, 1)
            if currCellFaces[faceI] == -1
                currCellFaces[faceI] = nFaces
                faces[nFaces] = cellFaces[cellI][faceI]

                nFaces += 1
            end
        end
    end
    #println(cellsAsFaces)

    resize!(faces, nFaces-1)

    #println(faces)
    return patchSizes, patchStarts, defaultPatchStart,
           faces, nFaces, cellsAsFaces

end

function calc_merge_info(blocks::Vector{Block})

    println("Creating block offsets")

    nBlocs = size(blocks, 1)
    blockOffsets = Vector{Int64}(nBlocks)
    
    nPoints = 0
    nCells = 0

    for blockI in 1:nBlocks
        blockOffsets[blockI] = nPoints
        nPoints += size(blocks[blockI].points, 1)
        nCells += size(blocks[blockI].cells, 1)
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

# Create mesh form the blocks
pointCellAddressing = point_cell_addressing(convert(Vector{Cell}, blocks), size(vertices, 1))
blockFaces = [cellfaces(convert(Cell, blocks[i])) for i in 1:nBlocks]
patchFaceCells = patch_face_cells(patchBlockFaces[1], blockFaces, pointCellAddressing)

blocksAsCells =Vector{Cell}(nBlocks)

for i in 1:nBlocks
    blocksAsCells[i] = convert(Cell, blocks[i])
end

patchSizes, patchStarts, defaultPatchStart, faces, nFaces, cellsAsFaces =
    create_topology(blocksAsCells,
                    patchBlockFaces,
                    patchNames,
                    size(vertices, 1))

nDefaultFaces = nFaces - defaultPatchStart

if nDefaultFaces > 0
    warn("Undefined block faces present in the mesh description")
end

owner = fill(-1, size(faces, 1))
neighbour = fill(-1, size(faces, 1))

markedFaces = fill(false, size(faces, 1))

nInternalFaces = 0

for cellI in 1:size(cellsAsFaces, 1)
    
    cellFaces = cellsAsFaces[cellI]

    for faceI in 1:size(cellFaces, 1)
        if cellFaces[faceI] < 1
            faceLabel = cellFaces[cellI]
            error("Illegal face label $faceLabel in cell $cellI")
        end

        if !markedFaces[cellFaces[faceI]]
            # First visit: owner
            owner[cellFaces[faceI]] = cellI
            markedFaces[cellFaces[faceI]] = true
        else
            # Second visit: neighbour
            neighbour[cellFaces[faceI]] = cellI
            nInternalFaces += 1
        end
    end
end

resize!(neighbour, nInternalFaces)

calc_merge_info(blocks)
end
