module blockMesh

include("MeshPrimitives.jl")
include("Edges.jl")
include("Blocks.jl")

using JSON
using .MeshPrimitives
using .Edges
using .Blocks


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


"""
    point_cell_addressing(cells::Vector{Cell}, nPoints::Int64))

Compute point to cell addressing.

Creates a `Vector` of `Vector`s  of size nPoints, each corresponding to a point
in the geometry. Each `Vector` contains the numbers of the cells that contain
current point.
"""
function point_cell_addressing(cells::Vector{Cell}, nPoints::Int64)

    # Number of cells each point is included in
    pointCellAddressing = Vector{Vector{Int64}}(nPoints)

    for i in 1:nPoints
        pointCellAddressing[i] = Vector{Int64}(0)
    end

    # For each cell
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

            # For each cell that this point is part of
            for cellI in 1:size(facePointCells, 1)

                cellIFaces = cellFaces[facePointCells[cellI]]

                # For each face of the cells check if it is the same as faceI
                for cellIFaceI in 1:size(cellIFaces, 1)

                    if samepoints(cellIFaces[cellIFaceI], faces[faceI])
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
                         boundaryPatchNames::Vector{String},
                         pointCellAddressing::Vector{Vector{Int64}},
                         nPoints)

    # Define a vector of cells defined as face index vectors
    cellsAsFaces = Vector{Vector{Int64}}(size(cells, 1))

    # Get the faces of each cell and the maximum number of faces
    cellFaces = Vector{Vector{Face}}(size(cells, 1))
    maxFaces = 0

    for i in 1:size(cellFaces,1)
        # get the faces of the cell
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

    found = false

    # Add the non-boundary faces of all cells to the face list
    # Form the cellsAsFaces array
    for cellI in 1:size(cells, 1)

        cellIFaces = cellFaces[cellI]

        # Record the neighbour cell
        neiCells = fill(-1, size(cellIFaces))

        # Record the face of neighbour cell
        faceOfNeiCell = fill(-1, size(cellIFaces))

        nNeighbours = 0

        # Identify the neighbours for each face of the cell
        for faceI in 1:size(cellIFaces, 1)

            # Skip faces that have already been matched
            if cellsAsFaces[cellI][faceI] > -1
                continue
            end

            found = false

            # the current face we are matching    
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
                # All cells should have at least one neighbour
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
                warn("Putting faces into default patch")
                currCellFaces[faceI] = nFaces
                faces[nFaces] = cellFaces[cellI][faceI]

                nFaces += 1
            end
        end
    end

    resize!(faces, nFaces-1)

    return patchSizes, patchStarts, defaultPatchStart,
           faces, nFaces-1, cellsAsFaces

end

function calc_merge_info(blocks::Vector{Block},
                         blockPoints::Vector{Point},
                         blockFaces::Vector{Face},
                         blockCellsAsFaces::Vector{Vector{Int64}},
                         faceOwnerBlocks::Vector{Int64},
                         faceNeighbourBlocks::Vector{Int64},
                         nInternalFaces::Int64)

    println("Creating block offsets")

    nBlocks = size(blocks, 1)
    blockOffsets = Vector{Int64}(nBlocks)

    nPoints = 0
    nCells = 0

    for blockI in 1:nBlocks
        blockOffsets[blockI] = nPoints
        nPoints += size(blocks[blockI].points, 1)
        nCells += size(blocks[blockI].cells, 1)
    end

    mergeList = fill(-1, nPoints)

    glueMergePairs = Vector{Vector{Vector{Int64}}}(size(blockFaces, 1))

    println("Creating merge list")

    # Go through all the faces in the topology
    for sI in 1:size(blockFaces, 1)

        # Grab the owner of the face
        blockPLabel = faceOwnerBlocks[sI]
        # Grab the points of the owner of the face
        blockPPoints = blocks[blockPLabel].points
        # Grab the faces of the onwer of the face
        blockPFaces = blockCellsAsFaces[blockPLabel]

        # Insure that one of the faces of the owner-block
        # is in fact sI
        foundFace = false

        blockPfaceLabel = 0
        for i in 1:size(blockPFaces, 1)
            if blockPFaces[i] == sI
                foundFace = true
                blockPfaceLabel = i
                break
            end
        end

        if !foundFace
            error("calc_merge_info(): cannot find merge face
                   for block $blockPLabel")
        end
    
        # faces on the block boundary (face)
        blockPfaceFaces = blocks[blockPLabel].boundaryFaces[blockPfaceLabel]

        glueMergePairs[sI] = [Vector{Int64}(1) 
                              for _ in 1:size(blockPfaceFaces, 1)]
        
        bbMin, bbMax = bounding_box(blocks[blockPLabel].points)
        mergeSqrDist = norm((bbMax-bbMin)*1e-5)^2

        sqrMergeTol = 1e+6

        # for all faces on the boundary
        for fI in 1:size(blockPfaceFaces, 1)
            
            # points of current face
            blockPfaceFacePoints = blockPfaceFaces[fI]

           
            # for all points of current face
            for pI in 1:size(blockPfaceFacePoints, 1)
                #compare to the other points on this face
                for pI2 in 1:size(blockPfaceFacePoints, 1)
                    if pI != pI2
                    
                        magSqrDist = norm(
                            blockPPoints[blockPfaceFacePoints[pI]] -
                            blockPPoints[blockPfaceFacePoints[pI2]])^2
                        if magSqrDist < mergeSqrDist
                            PpointLabel = 
                                blockPfaceFacePoints[pI] +
                                blockOffsets[blockPlabel];

                            PpointLabel2 =
                                blockPfaceFacePoints[pI2] +
                                blockOffsets_[blockPlabel];
                            
                            minPP2 = min(PpointLabel, PpointLabel2);

                            if mergeList[PpointLabel] != -1
                                minPP2 = min(minPP2, mergeList[PpointLabel]);
                            end

                            if (mergeList_[PpointLabel2] != -1)
                                minPP2 = min(minPP2, mergeList[PpointLabel2]);
                            end 
                            
                            mergeList_[PpointLabel2] = minPP2
                            mergeList[PpointLabel] = minPP2
                        else
                            sqrMergeTol = min(sqrMergeTol, magSqrDist)
                        end
                    end
                end
            end 
        end # for all faces on current blockface

        sqrMergeTol /= 10

        # if face is internal
        if sI <= nInternalFaces
            blockNLabel = faceNeighbourBlocks[sI];
            blockNPoints = blocks[blockNLabel].points;
            blockNFaces = blockCellsAsFaces[blockNLabel]; 
            
            blockNfaceLabel = 1
            foundFace = false
            for i in 1:size(blockNFaces, 1)
                if blockNFaces[i] == sI 
                   blockNfaceLabel = i    
                   foundFace = true     
                   break
               end
            end

            if !foundFace
                error("calc_merge_info(): cannot find merge face
                      for block $blockNLabel")
            end


            blockNfaceFaces =
                blocks[blockNLabel].boundaryFaces[blockNfaceLabel]

            #println(blockPfaceLabel)
            #println(blockNfaceLabel)
            #println(blockPfaceFaces)
            #println(blockNfaceFaces)

            if size(blockPfaceFaces, 1) != size(blockNfaceFaces, 1)     
                error("calc_merge_info(): inconsistent number of faces
                      between block par $blockPlabel $blockNlabel")
            end

            # n^2 point search over all points of all faces of
            # master block over all points of all faces of slave block
            for fI in 1:size(blockPfaceFaces, 1)
                blockPfaceFacePoints = blockPfaceFaces[fI]

                resize!(glueMergePairs[sI][fI], 
                        size(blockPfaceFacePoints, 1))
                #println("glueMergePairs $glueMergePairs")
                glueMergePairs[sI][fI] = 
                    fill(-1, size(blockPfaceFacePoints, 1))

                for pI in 1:size(blockPfaceFacePoints, 1)
                    
                    for fIN in 1:size(blockNfaceFaces, 1)
                        blockNfaceFacePoints =
                            blockNfaceFaces[fIN]

                        for pIN in 1:size(blockNfaceFacePoints, 1)         
                                #println(blockPPoints[blockPfaceFacePoints[pI]])
                                #println(blockNPoints[blockNfaceFacePoints[pIN]])
                                #println("")
                            if norm(
                                blockPPoints[blockPfaceFacePoints[pI]] -
                                blockNPoints[blockNfaceFacePoints[pIN]])^2 < sqrMergeTol
                                # Found a pair
                                #println("Found pair")

                                glueMergePairs[sI][fI][pI] =
                                    blockNfaceFacePoints[pIN]
                                PPointLabel =
                                    blockPfaceFacePoints[pI] + blockOffsets[blockPLabel]
                                    
                                NPointLabel =
                                    blockNfaceFacePoints[pIN] + blockOffsets[blockNLabel]
                                
                                minPN = min(PPointLabel, NPointLabel)    

                                if mergeList[PPointLabel] != -1
                                    minPN = min(minPN, mergeList[PPointLabel])
                                end

                                if mergeList[NPointLabel] != -1
                                    minPN = min(minPN, mergeList[NPointLabel])
                                end

                                mergeList[PPointLabel] = minPN
                                mergeList[NPointLabel] = minPN
                            end
                        end
                    end
                end # for point in current face of current blockface

                for pI in size(blockPfaceFacePoints, 1)
                    if glueMergePairs[sI][fI][pI] == -1
                        error("Inconsitent point location between block pair
                            $blockPLabel and $blockNLabel 
                            (probably due to inconsistent grading)")
                            
                    end
                end
            end # for face in current blockface
        end #if face is internal

     end # for all faces of blocks

     blockInternalFaces = blockFaces[1:nInternalFaces]
     changedPointMerge = true 

     nPasses = 0

     while changedPointMerge
         changedPointMerge = false
         nPasses += 1

         for sI in 1:size(blockInternalFaces, 1)
             blockPlabel = faceOwnerBlocks[sI]
             blockNlabel = faceNeighbourBlocks[sI]

             blockPfaces = blockCellsAsFaces[blockPlabel]
             blockNfaces = blockCellsAsFaces[blockNlabel]

             blockPfaceLabel = 0    
             for i in 1:size(blockPfaces, 1)
                 if blockFaces[blockPfaces[i]] ==
                    blockInternalFaces[sI]
                 
                    blockPfaceLabel = i
                    break
                end
             end

            blockNfaceLabel = 0    
            for i in 1:size(blockNfaces, 1)
                if blockFaces[blockNfaces[i]] ==
                   blockInternalFaces[sI]
                    blockNfaceLabel = i
                    break
                end
            end

            blockPfaceFaces = blocks[blockPlabel].boundaryFaces[blockPfaceLabel]

            for fI in size(blockPfaceFaces, 1)
            
                blockPfaceFacePoints = blockPfaceFaces[fI]

                
                for pI = 1:size(blockPfaceFacePoints, 1)
                    PpointLabel = 
                        blockPfaceFacePoints[pI] + blockOffsets[blockPlabel]
                        

                    NpointLabel = 
                        glueMergePairs[sI][fI][pI] + blockOffsets[blockNlabel]
                        

                    if mergeList[PpointLabel] != mergeList[NpointLabel]
                        changePointMerge = true
                        minPN = min(mergeList[PpointLabel], mergeList[NpointLabel])
                        mergeList[PpointLabel] = minPN
                        mergeList[NpointLabel] = minPN
                    end
                end
            end

         end # for each interal face

         if nPasses > 100
             error("calc_merge_info(): Point merging faild after max number
                   of passes.")
         end

     end # while

     for sI in 1:size(blockInternalFaces, 1)
        blockPlabel = faceOwnerBlocks[sI]
        blockNlabel = faceNeighbourBlocks[sI]

        blockPfaces = blockCellsAsFaces[blockPlabel]
        blockNfaces = blockCellsAsFaces[blockNlabel]

        blockPpoints = blocks[blockPlabel].points
        blockNpoints = blocks[blockNlabel].points

        foundFace = false
        blockPfaceLabel = 0

        for i in 1:size(blockPfaces, 1)
            if blockFaces[blockPfaces[i]] == blockInternalFaces[sI]
                blockPfaceLabel = i
                foundFace = true
            end
        end

        if !foundFace
            error("Cannot find merge face for block $blockPlabel")    
        end

        foundFace = false
        blockNfaceLabel = 0

        for i in 1:size(blockNfaces, 1)
            if blockFaces[blockNfaces[i]] == blockInternalFaces[sI]
                blockNfaceLabel = i
                foundFace = true
            end
        end

        if !foundFace
            error("Cannot find merge face for block $blockNlabel")    
        end

        blockPfaceFaces = blocks[blockPlabel].boundaryFaces[blockPfaceLabel]
        blockNfaceFaces = blocks[blockNlabel].boundaryFaces[blockNfaceLabel]


        for fI in 1:size(blockPfaceFaces, 1)
            blockPfaceFacePoints = blockPfaceFaces[fI]

            for pI in 1:size(blockPfaceFacePoints, 1)
                PpointLabel = blockPfaceFacePoints[pI] + 
                              blockOffsets[blockPlabel]

                if mergeList[PpointLabel] == -1              
                    error("Unable to merge point
                          $pI
                          of  face
                          $fI
                          of block
                          $blockPlabel")
                end
            end
        end

        for fIN in 1:size(blockNfaceFaces, 1)
            blockNfaceFacePoints = blockNfaceFaces[fIN]

            for pIN in 1:size(blockNfaceFacePoints, 1)
                NpointLabel = blockNfaceFacePoints[pIN] + 
                              blockOffsets[blockNlabel]

                if mergeList[NpointLabel] == -1              
                    error("Unable to merge point
                          $pIN
                          of  face
                          $fIN
                          of block
                          $blockNlabel")
                end
            end
        end
     end # for internal blockface

     #Sort merge list to return new point label (in new shorter list)
     newPointLabel = 1

     for pointLabel in 1:size(mergeList, 1)
        if mergeList[pointLabel] > pointLabel
            error("Merge list contains point index out of range")
        end

        if mergeList[pointLabel] == -1 || mergeList[pointLabel] == pointLabel
            mergeList[pointLabel] = newPointLabel
            newPointLabel += 1
        else
            mergeList[pointLabel] = mergeList[mergeList[pointLabel]]
        end
     end

     nPoints = newPointLabel - 1

     return nCells, nPoints, blockOffsets, mergeList

end

function create_points(blocks, blockOffsets, mergeList, nPoints)
    println("Creating global point list")
    points = Vector{Point}(nPoints)

    for blockI in 1:size(blocks, 1)
        blockPoints = blocks[blockI].points

        for blockPointI in 1:size(blockPoints, 1)
            points[mergeList[blockOffsets[blockI] + blockPointI]] = 
                blockPoints[blockPointI]
        end
    end

    return points
end

function create_cells(blocks, blockOffsets, mergeList, nCells)
    println("Creating global cell list")
    cells = Vector{Cell}(nCells)

    cellLabel = 1

    for blockI in 1:size(blocks, 1)
        blockCells = blocks[blockI].cells

        cellPoints = Vector{Int64}(0)
        for blockCellI in 1:size(blockCells,1 )

            resize!(cellPoints, size(blockCells[blockCellI], 1))

            for cellPointI in 1:size(cellPoints, 1)
                idx =  blockCells[blockCellI][cellPointI] +
                       blockOffsets[blockI]
                cellPoints[cellPointI] = mergeList[idx]    
            end

            cells[cellLabel] = cellPoints 
            cellLabel += 1
        end

    end
    return cells
end

function create_patches(blocks, blockOffsets, mergeList, blockFaces, patchTopologyFaces, faces, owner)

    println("Creating patches")

    patches = [Vector{Face}(0) for _ in 1:size(patchTopologyFaces, 1)]

    # compute the faces of each patch 
    for patchI in 1:size(patches, 1)
        
        # find the owners of the surfaces definin gthe patch
        blockLabels = Vector{Int64}(size(patchTopologyFaces[patchI], 1)) 

        for fI in 1:size(patchTopologyFaces[patchI], 1)

            foundOwner = false
            for faceI in 1:size(faces, 1)

                if samepoints(faces[faceI], patchTopologyFaces[patchI][fI])
                    foundOwner = true
                    blockLabels[fI] = owner[faceI]
                end
            end

            if foundOwner == false
                error("Could not find owner of face $fI of patch $patchI")
            end
        end

        nFaces = 0
        
        # compute the number of faces on the patch
        for patchTopologyFaceLabel in 1:size(patchTopologyFaces[patchI], 1)
            blockI = blockLabels[patchTopologyFaceLabel]

            blockFacesI = blockFaces[blockI]

            for blockFaceLabel in 1:size(blockFacesI, 1)
                if samepoints(blockFacesI[blockFaceLabel], 
                   patchTopologyFaces[patchI][patchTopologyFaceLabel]) 
                   
                    nFaces += size(blocks[blockI].boundaryFaces[blockFaceLabel], 1)
                end
            end
        end

        patchFaces = Vector{Face}(nFaces)
        faceLabel = 1

        quadFace = Vector{Int64}(4)

        for patchTopologyFaceLabel in 1:size(patchTopologyFaces[patchI], 1)

            blockI = blockLabels[patchTopologyFaceLabel]
            blockFacesI = blockFaces[blockI]

            for blockFaceLabel in 1:size(blockFacesI, 1)

                if samepoints(blockFacesI[blockFaceLabel], 
                   patchTopologyFaces[patchI][patchTopologyFaceLabel]) 

                    blockPatchFaces = blocks[blockI].boundaryFaces[blockFaceLabel]

                    for fI in 1:size(blockPatchFaces, 1)
                        idx = blockPatchFaces[fI][1] +
                              blockOffsets[blockI]
                        quadFace[1] = mergeList[idx]      
                              
                        nUnique = 2

                        for facePointLabel in 2:4

                            idx = blockPatchFaces[fI][facePointLabel] +
                                  blockOffsets[blockI]
                            quadFace[nUnique] = mergeList[idx]

                            if quadFace[nUnique] != quadFace[nUnique - 1]
                                nUnique += 1
                            end

                        end

                        if quadFace[nUnique - 1] == quadFace[1]
                            nUnique -= 1
                        end

                        # re-adjust due to fiddling to get 1-based array  to
                        # work
                        nUnique -= 1

                        if nUnique == 4
                            patchFaces[faceLabel] = quadFace
                            faceLabel += 1
                        elseif nUnuque == 3
                            warn("Boundary face does not have 4 unique points")    
                            patchFaces[faceLabel] = quadFace[1:3]
                            faceLabel += 1

                        end
                    end
                end #if same face
            end
        end

        resize!(patchFaces, faceLabel - 1)
        patches[patchI] = patchFaces

    end # for patchI


    return patches
end

function init_mesh(faces, cellsAsFaces)
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

    return owner, neighbour, nInternalFaces
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
        println("    Block number ", blockI)
        blocks[blockI] = Block()
        blocks[blockI].nCells = dict["blocks"][blockI]["number of cells"]

        # read vertex numbers defining the block
        blocks[blockI].vertexLabels = dict["blocks"][blockI]["vertices"] + 1

        # the coordinates of all vertices defining the block
        # corresponding to the vertex numbers defining the block
        for j in 1:8
            blocks[blockI].vertices[j] = vertices[blocks[blockI].vertexLabels[j]]
        end

        println("        Creating edge-points")
        make_block_edges!(blocks[blockI])

        println("        Creating points")
        create_points!(blocks[blockI])

        println("        Creating cells")
        create_cells!(blocks[blockI])

        println("        Creating boundary faces")
        create_boundary_faces!(blocks[blockI])
        println("        Done")
    end

    return blocks
end

function main()
    dictPath = joinpath("tests", "dict.json")

    # Parse the dictionary
    dict = JSON.parsefile(dictPath)

    # Get the number of blocks
    nBlocks = size(dict["blocks"], 1)
    nVertices = size(dict["vertices"], 1)

    #vertices defining the mesh as defined in the dict
    vertices = Vector{Point}(nVertices)
    vertices[:] = dict["vertices"][:]

    patchNames = String[]
    # The faces of the patches defined from vertices in the dict
    patchSurfaces = Vector{Vector{Face}}(0)

    patchNames, patchSurfaces =  read_boundary!(dict)
    check_patch_vertex_labels(patchNames, patchSurfaces, vertices)

    println("Creating blocks...")
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


    println("Writing")
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


end
end
