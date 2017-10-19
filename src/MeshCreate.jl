
export point_cell_addressing, patch_face_cells, create_points, create_cells,
create_patches, create_topology, calc_merge_info, create_owner_neighbour 


"""
    point_cell_addressing(cells::Vector{Cell}, nPoints::Label))

Compute point to cell addressing.

Creates a `Vector` of `Vector`s  of size nPoints, each corresponding to a point
in the geometry. Each `Vector` contains the numbers of the cells that contain
current point.
"""
function point_cell_addressing(
    cells::Vector{Cell{Label}},
    nPoints::Label
) where {Label <: Union{Int32, Int64}}
    @inbounds begin

    # Number of cells each point is included in
    pointCellAddressing = Vector{Vector{Label}}(nPoints)

    for i in 1:nPoints
        pointCellAddressing[i] = Vector{Label}(0)
    end

    # For each cell
    for cellI in eachindex(cells)

        # For each point in the cell
        for pointI in eachindex(cells[cellI])

            currPoint = cells[cellI][pointI]
            # Add current cell to the addressing list
            push!(pointCellAddressing[currPoint], cellI)
        end
    end

    end #inbounds
    return pointCellAddressing
end

function patch_face_cells(
    faces::Vector{Face{Label}},
    cellFaces::Vector{Vector{Face{Label}}},
    pointCellAddressing::Vector{Vector{Label}}
) where {Label <: Union{Int32, Int64}}

    @inbounds begin

    faceCells = Vector{Label}(length(faces))

    # For each face determine the cell
    for faceI in eachindex(faces)

        found = false
        currFace = faces[faceI]

        # For each point of the face
        for facePointI in eachindex(currFace)
            facePointCells = pointCellAddressing[currFace[facePointI]]

            # For each cell that this point is part of
            for cellI in eachindex(facePointCells)

                cellIFaces = cellFaces[facePointCells[cellI]]

                # For each face of the cells check if it is the same as faceI
                for cellIFaceI in eachindex(cellIFaces)

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

    end #inbounds

    return faceCells
end


function create_topology(
    cells::Vector{Cell{Label}},
    boundaryFaces::Vector{Vector{Face{Label}}},
    boundaryPatchNames::Vector{String},
    pointCellAddressing::Vector{Vector{Label}},
    nPoints::Label
) where {Label<:Union{Int32, Int64}}
    @inbounds begin

    # Define a vector of cells defined as face index vectors
    cellsAsFaces = Vector{Vector{Label}}(length(cells))

    # Get the faces of each cell and the maximum number of faces
    cellFaces = Vector{Vector{Face{Label}}}(length(cells))
    maxFaces::Label = 0

    for i in eachindex(cellFaces)
        # get the faces of the cell
        cellFaces[i] = cellfaces(cells[i])
        maxFaces += length(cellFaces[i])

        # Set the associated face label to -1, to mark as undefined
        cellsAsFaces[i] = fill(-1, length(cellFaces[i]))
    end

    # Declare the array of faces
    faces = Vector{Face{Label}}(maxFaces)

    nFaces::Label = 1

    # Get point to cell addressing
    pointCellAddr = point_cell_addressing(cells, nPoints)

    found = false

    # Add the non-boundary faces of all cells to the face list
    # Form the cellsAsFaces array
    for cellI in eachindex(cells)

        cellIFaces = cellFaces[cellI]

        # Record the neighbour cell
        neiCells = fill(Label(-1), length(cellIFaces))

        # Record the face of neighbour cell
        faceOfNeiCell = fill(Label(-1), length(cellIFaces))

        nNeighbours::Label = 0

        # Identify the neighbours for each face of the cell
        for faceI in eachindex(cellIFaces)

            # Skip faces that have already been matched
            if cellsAsFaces[cellI][faceI] > -1
                continue
            end

            found = false

            # the current face we are matching    
            currFace = cellIFaces[faceI]

            for pointI in eachindex(currFace)

                # Get list of cells sharing this point
                currNeighbours = pointCellAddr[currFace[pointI]]

                for neiI in eachindex(currNeighbours)

                    currNei = currNeighbours[neiI]

                    # Reject neighbours with lower index
                    if currNei > cellI

                        # The list of faces to search through
                        neiFaces = cellFaces[currNei]

                        for neiFaceI in eachindex(neiFaces)

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
            nextNei::Label = -1
            minNei::Label = length(cellsAsFaces)

            for ncI in eachindex(neiCells)
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

    patchSizes = Vector{Label}(length(boundaryFaces))
    patchStarts = Vector{Label}(length(boundaryFaces))

    for patchI in eachindex(boundaryFaces)

        patchFaces = boundaryFaces[patchI]
        patchName = boundaryPatchNames[patchI]

        currPatchFaceCells = patch_face_cells(patchFaces, cellFaces,
                                              pointCellAddressing)
        currPatchStart::Label = nFaces


        for faceI in eachindex(patchFaces)

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

    defaultPatchStart::Label = nFaces

    # Take care of "non-existing faces", put them into the default patch
    for cellI in eachindex(cellsAsFaces)

        currCellFaces = cellsAsFaces[cellI]

        for faceI in eachindex(currCellFaces)
            if currCellFaces[faceI] == -1
                warn("Putting faces into default patch")
                currCellFaces[faceI] = nFaces
                faces[nFaces] = cellFaces[cellI][faceI]

                nFaces += 1
            end
        end
    end

    nFaces -= 1
    resize!(faces, nFaces)

    end #inbounds 

    return patchSizes, patchStarts, defaultPatchStart,
           faces, nFaces, cellsAsFaces

end


function calc_merge_info(
    blocks::Vector{Block{Label}},
    blockPoints::Vector{Point},
    blockFaces::Vector{Face{Label}},
    blockCellsAsFaces::Vector{Vector{Label}},
    faceOwnerBlocks::Vector{Label},
    faceNeighbourBlocks::Vector{Label},
    nInternalFaces::Label
) where {Label <: Union{Int32, Int64}}
    @inbounds begin

    nBlocks::Label = length(blocks)
    blockOffsets = Vector{Label}(nBlocks)

    nPoints::Label = 0
    nCells::Label = 0

    for blockI in 1:nBlocks
        blockOffsets[blockI] = nPoints
        nPoints += length(blocks[blockI].points)
        nCells += length(blocks[blockI].cells)
    end

    mergeList = fill(Label(-1), nPoints)

    glueMergePairs = Vector{Vector{Vector{Label}}}(length(blockFaces))


    # Go through all the faces in the topology
    for sI in eachindex(blockFaces)

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
    
        # faces on the block surface
        blockPfaceFaces = blocks[blockPLabel].boundaryFaces[blockPfaceLabel]

        glueMergePairs[sI] = [Vector{Label}(1) 
                              for _ in 1:size(blockPfaceFaces, 1)]
        
        bbMin, bbMax = bounding_box(blocks[blockPLabel].points)
        mergeSqrDist = norm((bbMax-bbMin)*1e-5)^2

        sqrMergeTol = 1e+6

        # for all faces on the sirface
        for fI in eachindex(blockPfaceFaces)
            
            # points of current face
            blockPfaceFacePoints = blockPfaceFaces[fI]

           
            # for all points of current face
            for pI in eachindex(blockPfaceFacePoints)
                #compare to the other points on this face
                for pI2 in eachindex(blockPfaceFacePoints)
                    if pI != pI2
                    
                        magSqrDist = norm(
                            blockPPoints[blockPfaceFacePoints[pI]] -
                            blockPPoints[blockPfaceFacePoints[pI2]])^2
                        if magSqrDist < mergeSqrDist
                            PpointLabel = 
                                blockPfaceFacePoints[pI] +
                                blockOffsets[blockPLabel];

                            PpointLabel2 =
                                blockPfaceFacePoints[pI2] +
                                blockOffsets[blockPLabel];
                            
                            minPP2 = min(PpointLabel, PpointLabel2);

                            if mergeList[PpointLabel] != -1
                                minPP2 = min(minPP2, mergeList[PpointLabel]);
                            end

                            if (mergeList[PpointLabel2] != -1)
                                minPP2 = min(minPP2, mergeList[PpointLabel2]);
                            end 
                            
                            mergeList[PpointLabel2] = minPP2
                            mergeList[PpointLabel] = minPP2
                        else
                            sqrMergeTol = min(sqrMergeTol, magSqrDist)
                        end
                    end
                end
            end 
        end # for all faces on current blockface

        sqrMergeTol /= 10

        # if suraface is internal
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

                glueMergePairs[sI][fI] = 
                    fill(-1, size(blockPfaceFacePoints, 1))

                for pI in 1:size(blockPfaceFacePoints, 1)
                    
                    for fIN in 1:size(blockNfaceFaces, 1)
                        blockNfaceFacePoints =
                            blockNfaceFaces[fIN]

                        for pIN in 1:size(blockNfaceFacePoints, 1)         

                            p1 = blockPPoints[blockPfaceFacePoints[pI]] 
                            p2 = blockNPoints[blockNfaceFacePoints[pIN]]

                            diff = p1 - p2
                            sqrNorm = diff[1]*diff[1] +
                                      diff[2]*diff[2] +
                                      diff[3]*diff[3]

                            if sqrNorm < sqrMergeTol
                                # Found a pair

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
                            $(blockPLabel) and $(blockNLabel)
                            (probably due to inconsistent grading)")
                            
                    end
                end
            end # for face in current blockface
        end #if surface is internal

     end # for all surfaces of blocks

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
     
     end # inbounds
     return nCells, nPoints, blockOffsets, mergeList
     
end

function create_points(
    blocks,
    blockOffsets,
    mergeList,
    nPoints
)

    @inbounds begin

    points = Vector{Point}(nPoints)

    for blockI in 1:size(blocks, 1)
        blockPoints = blocks[blockI].points

        for blockPointI in 1:size(blockPoints, 1)
            points[mergeList[blockOffsets[blockI] + blockPointI]] = 
                blockPoints[blockPointI]
        end
    end

    end #indbounds
    return points
end

function create_cells(
    blocks::Vector{Block{Label}},
    blockOffsets::Vector{Label},
    mergeList::Vector{Label},
    nCells::Label
) where {Label <: Union{Int32, Int64}}

    @inbounds begin

    cells = Vector{Cell{Label}}(nCells)

    cellLabel::Label = 1

    for blockI in 1:size(blocks, 1)
        blockCells = blocks[blockI].cells

        cellPoints = Vector{Label}(0)
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

    end #inbounds
    return cells
end

function create_patches(
    blocks::Vector{Block{Label}},
    blockOffsets::Vector{Label},
    mergeList::Vector{Label},
    blockFaces::Vector{Vector{Face{Label}}},
    patchTopologyFaces::Vector{Vector{Face{Label}}},
    faces::Vector{Face{Label}},
    owner::Vector{Label}
) where {Label <: Union{Int32, Int64}}
    @inbounds begin

    patches = [Vector{Face{Label}}(0) for _ in 1:size(patchTopologyFaces, 1)]

    # compute the faces of each patch 
    for patchI in eachindex(patches)
        
        # find the owners of the surfaces definin gthe patch
        blockLabels = Vector{Label}(size(patchTopologyFaces[patchI], 1)) 

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

        nFaces::Label = 0
        
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
        faceLabel::Label = 1

        quadFace = Vector{Label}(4)

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
                              
                        nUnique::Label = 2

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
                        elseif nUnique == 3
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


    end #inbounds
    return patches
end

function create_owner_neighbour(
    nFaces::Label,
    cellsAsFaces::Vector{Vector{Label}}
) where {Label <: Union{Int32, Int64}}
   
    @inbounds begin

    owner = fill(Label(-1), nFaces)
    neighbour = fill(Label(-1), nFaces)

    markedFaces = fill(false, nFaces)

    nInternalFaces::Label = 0

    for cellI in eachindex(cellsAsFaces)

        cellFaces = cellsAsFaces[cellI]

        for faceI in eachindex(cellFaces)
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

    end #inbounds
    return owner, neighbour
end
