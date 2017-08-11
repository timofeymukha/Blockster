
using StaticArrays: SVector

#export Point, Face, Cell, cellfaces, samepoints


" Type for a Point in 3d space. An alias for a vector of Float64 of size 3"
const Point = SVector{3, Float64}

"""
Type for a quad-face. Holds a list of Point-labels. An alias for a vector of 
Int64 of size 4.
"""
const Face = SVector{4, Int64 }

"""
Type for a hex-cell. Holds a list of Point-labels. An alias for a vector of 
Int64 of size 8.
"""
const Cell = SVector{8, Int64 }

"""
Return the faces of a cell.

**Parameters:**

    1. cell::Cell
    The cell to return the faces for.

**Returns:**
    
    1. Vector{face}
    A vector of faces.
"""
function cellfaces(c::Cell)
    faces = Vector{Face}(6)

    faces[1] = Face([c[4], c[1], c[5], c[8]])
    faces[2] = Face([c[3], c[7], c[6], c[2]])
    faces[3] = Face([c[6], c[5], c[1], c[2]])
    faces[4] = Face([c[7], c[3], c[4], c[8]])
    faces[5] = Face([c[1], c[2], c[3], c[4]])
    faces[6] = Face([c[5], c[6], c[7], c[8]])

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

        # Count the occurences of this Point in faceA
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

#end
