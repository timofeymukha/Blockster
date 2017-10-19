
export bounding_box, cellfaces, samepoints, check_patch_vertex_labels

" Type for a Point in 3d space. An alias for a vector of Float64 of size 3"
const Point = SVector{3, Float64}

"""
Type for a quad-face. Holds a list of Point-labels. An alias for a vector of 
Label of size 4.
"""
const Face{T} = SVector{4, T}


"""
Type for a hex-cell. Holds a list of Point-labels. An alias for a vector of 
Label of size 8.
"""
const Cell{T} = SVector{8, T}

"""
Return the cell as a list of faces.

**Parameters:**

    1. cell::Cell
    The cell to return as list of faces.

**Returns:**
    
    1. Vector{face}
    A vector of faces.
"""
function cellfaces(c::Cell{Label}) where {Label <: Integer}
    @inbounds begin
    faces = Vector{Face{Label}}(6)

    faces[1] = Face{Label}([c[4], c[1], c[5], c[8]])
    faces[2] = Face{Label}([c[3], c[7], c[6], c[2]])
    faces[3] = Face{Label}([c[6], c[5], c[1], c[2]])
    faces[4] = Face{Label}([c[7], c[3], c[4], c[8]])
    faces[5] = Face{Label}([c[1], c[2], c[3], c[4]])
    faces[6] = Face{Label}([c[5], c[6], c[7], c[8]])

    end #inbounds
    return faces
end

"""
Test whether two faces consist of the same Points.

**Parameters:**

    1. faceA::Face
    The first face.

    2. faceB::Face
    The second face.

**Returns:**

    1. Bool
    True if the two faces consist of the same vertices.
"""
function samepoints(faceA::Face, faceB::Face)

    for PointI in 1:4
        # Count the occurences of this point in faceA
        faceAOcc = 0
        for i in 1:4
            if faceA[PointI] == faceA[i]
                faceAOcc += 1
            end
        end

        # Count the occurences of this Point in faceB
        faceBOcc = 0
        for i in 1:4
            if faceA[PointI] == faceB[i]
                faceBOcc += 1
            end
        end

        if faceAOcc != faceBOcc
            return false
        end
    end

    return true
end

"""
    bounding_box(points)

    Compute the bounding box of a set of points.
"""
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
