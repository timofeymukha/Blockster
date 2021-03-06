
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
) where {Label <: Union{Int32,Int64}}
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
) where {Label <: Union{Int32,Int64}}

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
) where {Label<:Union{Int32,Int64}}
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
                currNeighbours = pointCellAddressing[currFace[pointI]]

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
           faces, cellsAsFaces

end


function calc_merge_info(
    blocks::Vector{Block{Label}},
    surfaces::Vector{Face{Label}},
    blocksAsSurfaces::Vector{Vector{Label}},
    ownerBlocks::Vector{Label},
    neighbourBlocks::Vector{Label},
    nInternalSurfaces::Label
) where {Label <: Union{Int32,Int64}}
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

    glueMergePairs = Vector{Vector{Vector{Label}}}(length(surfaces))


    # Go through all the surfaces in the topology
    for sI in eachindex(surfaces)

        # Grab the owner of the surface
        blockPLabel = ownerBlocks[sI]
        # Grab the points of the owner of the surface
        blockPPoints = blocks[blockPLabel].points
        # Grab the surfaces of the onwer of the surface
        blockPSurfaces = blocksAsSurfaces[blockPLabel]

        # Insure that one of the suraces of the owner-block
        # is in fact sI
        foundSurface = false

        blockPSurfaceLabel = 0
        for i in 1:size(blockPSurfaces, 1)
            if blockPSurfaces[i] == sI
                foundSurface = true
                blockPSurfaceLabel = i
                break
            end
        end

        if !foundSurface
            error("calc_merge_info(): cannot find merge face
                   for block $blockPLabel")
        end
    
        # faces on the block surface
        blockPSurfaceFaces = blocks[blockPLabel].boundaryFaces[blockPSurfaceLabel]

        glueMergePairs[sI] = [Vector{Label}(1) 
                              for _ in 1:size(blockPSurfaceFaces, 1)]
        
        bbMin, bbMax = bounding_box(blocks[blockPLabel].points)
        mergeSqrDist = norm((bbMax - bbMin)*1e-5)^2

        sqrMergeTol = 1e+6

        # for all faces on the surface
        for fI in eachindex(blockPSurfaceFaces)
            
            # points of current face
            blockPSurfaceFacePoints = blockPSurfaceFaces[fI]

           
            # for all points of current face
            for pI in eachindex(blockPSurfaceFacePoints)
                #compare to the other points on this face
                for pI2 in eachindex(blockPSurfaceFacePoints)
                    if pI != pI2
                    
                        magSqrDist = norm(
                            blockPPoints[blockPSurfaceFacePoints[pI]] -
                            blockPPoints[blockPSurfaceFacePoints[pI2]])^2

                        # Points are so close, we need to merge them
                        if magSqrDist < mergeSqrDist
                            PpointLabel = 
                                blockPSurfaceFacePoints[pI] +
                                blockOffsets[blockPLabel];

                            PpointLabel2 =
                                blockPSurfaceFacePoints[pI2] +
                                blockOffsets[blockPLabel];
                            
                            minPointLabel = min(PpointLabel, PpointLabel2);

                            if mergeList[PpointLabel] != -1
                                minPointLabel = min(minPP2, mergeList[PpointLabel]);
                            end

                            if (mergeList[PpointLabel2] != -1)
                                minPointLabel = min(minPP2, mergeList[PpointLabel2]);
                            end 
                            
                            mergeList[PpointLabel2] = minPointLabel
                            mergeList[PpointLabel] = minPointLabel
                        else
                            sqrMergeTol = min(sqrMergeTol, magSqrDist)
                        end
                    end
                end
            end 
        end # for all faces on current surface 

        sqrMergeTol /= 10

        # if suraface is internal
        if sI <= nInternalSurfaces
            blockNLabel = neighbourBlocks[sI];
            blockNPoints = blocks[blockNLabel].points;
            blockNSurfaces = blocksAsSurfaces[blockNLabel]; 
            
            blockNSurfaceLabel = 1
            foundSurface = false
            for i in 1:size(blockNSurfaces, 1)
                if blockNSurfaces[i] == sI 
                   blockNSurfaceLabel = i    
                   foundSurface = true     
                   break
               end
            end

            if !foundSurface
                error("calc_merge_info(): cannot find merge face
                      for block $blockNLabel")
            end


            blockNSurfaceFaces =
                blocks[blockNLabel].boundaryFaces[blockNSurfaceLabel]

            if length(blockPSurfaceFaces) != length(blockNSurfaceFaces)     
                error("calc_merge_info(): inconsistent number of faces
                      between block pair $blockPlabel $blockNlabel")
            end

            # n^2 point search over all points of all faces of
            # master block over all points of all faces of slave block

            # for each face on master block
            for fI in eachindex(blockPSurfaceFaces)
                blockPSurfaceFacePoints = blockPSurfaceFaces[fI]

                resize!(glueMergePairs[sI][fI], 
                        size(blockPSurfaceFacePoints, 1))

                glueMergePairs[sI][fI] = 
                    fill(Label(-1), length(blockPSurfaceFacePoints))

                # for each point on the face of the master block
                for pI in eachindex(blockPSurfaceFacePoints)
                    
                    # for each face of the slave block
                    for fIN in eachindex(blockNSurfaceFaces)
                        blockNSurfaceFacePoints =
                            blockNSurfaceFaces[fIN]

                        # for each point on the face of the slave block
                        for pIN in eachindex(blockNSurfaceFacePoints)         

                            p1 = blockPPoints[blockPSurfaceFacePoints[pI]] 
                            p2 = blockNPoints[blockNSurfaceFacePoints[pIN]]

                            diff = p1 - p2
                            sqrNorm = diff[1]*diff[1] +
                                      diff[2]*diff[2] +
                                      diff[3]*diff[3]

                            if sqrNorm < sqrMergeTol
                                # Found a pair

                                glueMergePairs[sI][fI][pI] =
                                    blockNSurfaceFacePoints[pIN]

                                PPointLabel =
                                    blockPSurfaceFacePoints[pI] + blockOffsets[blockPLabel]
                                    
                                NPointLabel =
                                    blockNSurfaceFacePoints[pIN] + blockOffsets[blockNLabel]
                                
                                minPointLabel = min(PPointLabel, NPointLabel)    

                                if mergeList[PPointLabel] != -1
                                    minPointLabel = 
                                        min(minPointLabel, mergeList[PPointLabel])
                                end

                                if mergeList[NPointLabel] != -1
                                    minPointLabel = 
                                        min(minPointLabel, mergeList[NPointLabel])
                                end

                                mergeList[PPointLabel] = minPointLabel
                                mergeList[NPointLabel] = minPointLabel
                            end
                        end
                    end
                end # for point in current face of current blockface

                for pI in eachindex(blockPSurfaceFacePoints)
                    if glueMergePairs[sI][fI][pI] == -1
                        error("Inconsitent point location between block pair
                            $(blockPLabel) and $(blockNLabel)
                            (probably due to inconsistent grading)")
                            
                    end
                end
            end # for face in current block surface
        end #if surface is internal
    end # for all surfaces of blocks

    blockInternalSurfaces = surfaces[1:nInternalSurfaces]
    changedPointMerge = true 

    nPasses = 0

    while changedPointMerge
        changedPointMerge = false
        nPasses += 1

        for sI in eachindex(blockInternalSurfaces)
            blockPlabel = ownerBlocks[sI]
            blockNlabel = neighbourBlocks[sI]

            blockPSurfaces = blocksAsSurfaces[blockPlabel]
            blockNSurfaces = blocksAsSurfaces[blockNlabel]

            blockPSurfaceLabel = 0    
            for i in 1:size(blockPSurfaces, 1)
                if surfaces[blockPSurfaces[i]] ==
                   blockInternalSurfaces[sI]
                
                   blockPSurfaceLabel = i
                   break
               end
            end

            blockNSurfaceLabel = 0    
            for i in 1:size(blockNSurfaces, 1)
                if surfaces[blockNSurfaces[i]] ==
                   blockInternalSurfaces[sI]
                    blockNSurfaceLabel = i
                    break
                end
            end

            blockPSurfaceFaces = blocks[blockPlabel].boundaryFaces[blockPSurfaceLabel]

            for fI in eachindex(blockPSurfaceFaces)
            
                blockPSurfaceFacePoints = blockPSurfaceFaces[fI]

                
                for pI in eachindex(blockPSurfaceFacePoints)
                    pPointLabel = 
                        blockPSurfaceFacePoints[pI] + blockOffsets[blockPlabel]
                        

                    nPointLabel = 
                        glueMergePairs[sI][fI][pI] + blockOffsets[blockNlabel]
                        

                    if mergeList[pPointLabel] != mergeList[nPointLabel]
                        changePointMerge = true
                        minPointLabel = min(mergeList[pPointLabel], mergeList[nPointLabel])
                        mergeList[pPointLabel] = minPointLabel
                        mergeList[nPointLabel] = minPointLabel
                    end
                end
            end

        end # for each interal face

        if nPasses > 100
            error("calc_merge_info(): Point merging failed after max number
                  of passes.")
        end

     end # while

    # Check all points have been merged
    for sI in eachindex(blockInternalSurfaces)
        blockPlabel = ownerBlocks[sI]
        blockNlabel = neighbourBlocks[sI]

        blockPSurfaces = blocksAsSurfaces[blockPlabel]
        blockNSurfaces = blocksAsSurfaces[blockNlabel]

        blockPpoints = blocks[blockPlabel].points
        blockNpoints = blocks[blockNlabel].points

        foundSurface = false
        blockPfaceLabel = 0

        for i in 1:size(blockPSurfaces, 1)
            if surfaces[blockPSurfaces[i]] == blockInternalSurfaces[sI]
                blockPfaceLabel = i
                foundSurface = true
            end
        end

        if !foundSurface
            error("Cannot find merge face for block $blockPlabel")    
        end

        foundSurface = false
        blockNSurfaceLabel = 0

        for i in 1:size(blockNSurfaces, 1)
            if surfaces[blockNSurfaces[i]] == blockInternalSurfaces[sI]
                blockNSurfaceLabel = i
                foundSurface = true
            end
        end

        if !foundSurface
            error("Cannot find merge face for block $blockNlabel")    
        end

        blockPSurfaceFaces = blocks[blockPlabel].boundaryFaces[blockPfaceLabel]
        blockNSurfaceFaces = blocks[blockNlabel].boundaryFaces[blockNSurfaceLabel]


        for fI in eachindex(blockPSurfaceFaces)
            blockPSurfaceFacePoints = blockPSurfaceFaces[fI]

            for pI in 1:size(blockPSurfaceFacePoints, 1)
                PpointLabel = blockPSurfaceFacePoints[pI] + 
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

        for fIN in eachindex(blockNSurfaceFaces)
            blockNSurfaceFacePoints = blockNSurfaceFaces[fIN]

            for pIN in eachindex(blockNSurfaceFacePoints)
                NpointLabel = blockNSurfaceFacePoints[pIN] + 
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
     end # for internal block surface

    #Sort merge list to return new point label (in new shorter list)
    newPointLabel = 1

    for pointLabel in eachindex(mergeList)
       if mergeList[pointLabel] > pointLabel
           error("Merge list contains point index out of range")
       end

       # if point is internal or had the least index of the two merged points
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

    for blockI in eachindex(blocks)
        blockPoints = blocks[blockI].points

        # Put each point according to its index in the merge list
        for blockPointI in eachindex(blockPoints)
            points[mergeList[blockOffsets[blockI] + blockPointI]] = 
                blockPoints[blockPointI]
        end
    end

    end #inbounds
    return points
end

function create_cells(
    blocks::Vector{Block{Label}},
    blockOffsets::Vector{Label},
    mergeList::Vector{Label},
    nCells::Label
) where {Label <: Union{Int32,Int64}}

    @inbounds begin

    cells = Vector{Cell{Label}}(nCells)

    cellLabel::Label = 1

    for blockI in eachindex(blocks)
        blockCells = blocks[blockI].cells

        for blockCellI in eachindex(blockCells)
            cellPoints = Vector{Label}(length(blockCells[blockCellI]))

            for cellPointI in eachindex(cellPoints)
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
    blocksAsSurfaces::Vector{Vector{Face{Label}}},
    patchSurfaces::Vector{Vector{Face{Label}}},
    surfaces::Vector{Face{Label}},
    ownerBlocks::Vector{Label}
) where {Label <: Union{Int32,Int64}}
    @inbounds begin

    patches = [Vector{Face{Label}}(0) for _ in 1:size(patchSurfaces, 1)]

    # compute the faces of each patch 
    for patchI in eachindex(patches)
        
        # find the owners of the surfaces defining this patch
        blockLabels = Vector{Label}(length(patchSurfaces[patchI])) 

        for patchSurfaceI in eachindex(patchSurfaces[patchI])

            foundOwner = false
            for surfaceI in eachindex(surfaces)

                if samepoints(surfaces[surfaceI], patchSurfaces[patchI][patchSurfaceI])
                    foundOwner = true
                    blockLabels[patchSurfaceI] = ownerBlocks[surfaceI]
                end
            end

            if foundOwner == false
                error("Could not find owner of facesIfI of patch $patchI")
            end
        end

        nFaces::Label = 0
        
        # compute the number of faces on the patch
        for patchSurfaceI in eachindex(patchSurfaces[patchI])

            # Get the label of the owner block
            blockI = blockLabels[patchSurfaceI]

            # The list of surfaces for the owner block
            blockAsSurfacesI = blocksAsSurfaces[blockI]

            # Among the list find the one coinciding with the
            # patch surface. Add the amount of boundary faces.
            for blockSurfaceI in eachindex(blockAsSurfacesI)
                if samepoints(blockAsSurfacesI[blockSurfaceI], 
                              patchSurfaces[patchI][patchSurfaceI]) 
                   
                    nFaces += size(blocks[blockI].boundaryFaces[blockSurfaceI], 1)
                end
            end
        end

        patchFaces = Vector{Face}(nFaces)
        faceLabel::Label = 1

        quadFace = Vector{Label}(4)

        for patchSurfaceFaceI in eachindex(patchSurfaces[patchI])

            blockI = blockLabels[patchSurfaceFaceI]
            blockAsSurfacesI = blocksAsSurfaces[blockI]

            for surfaceI in eachindex(blockAsSurfacesI)

                if samepoints(blockAsSurfacesI[surfaceI], 
                              patchSurfaces[patchI][patchSurfaceFaceI]) 

                    blockBoundaryFaces = blocks[blockI].boundaryFaces[surfaceI]

                    for fI in eachindex(blockBoundaryFaces)
                        idx = blockBoundaryFaces[fI][1] +
                              blockOffsets[blockI]
                        quadFace[1] = mergeList[idx]      
                              
                        nUnique::Label = 2

                        for facePointLabel in 2:4

                            idx = blockBoundaryFaces[fI][facePointLabel] +
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
                            error("Boundary face does not have 4 unique points")    
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

"""
    create_owner_neighbour(nFaces, cellsAsFaces)

Create owner and neigbour arrays.

Goes through the faces of each cell. The cell that contains a faces first
becomes its owner. The second cell to contain the same face becomes the 
neighbour. Note that the number of neighbours is equal to the number of 
internal faces.
"""
function create_owner_neighbour(
    nFaces::Label,
    cellsAsFaces::Vector{Vector{Label}}
) where {Label <: Union{Int32,Int64}}
   
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
