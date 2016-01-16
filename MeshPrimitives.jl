module MeshPrimitives

using FixedSizeArrays

export point, face, cell


" Type for a point in 3d space. An alias for a vector of Float64 of size 3"
typealias point Vec{3, Float64}

"""
Type for a quad-face. Holds a list of point-labels. An alias for a vector of 
Int64 of size 4.
"""
typealias face Vec{4, Int64 }

"""
Type for a hex-cell. Holds a list of point-labels. An alias for a vector of 
Int64 of size 8.
"""
typealias cell Vec{8, Int64 }

"""
Test whether two faces consist of the same points.

**Parameters:**

    1. faceA::face
    The first face.

    2. faceB::face
    The second face.

**Returns:**

    1. Bool
    True if the two faces consist of the same vertices.
"""
function samePoints(faceA::face, faceB::face)

    for pointI in 1:4
        # Count the occurences of this point in faceA
        faceAOcc = 0
        for i in 1:4
            if faceA[pointI] == faceA[i]
                faceAOcc += 1
            end
        end

        # Count the occurences of this point in faceA
        faceBOcc = 0
        for i in 1:4
            if faceB[pointI] == faceB[i]
                faceBOcc += 1
            end
        end

        if faceAOcc != faceBOcc
            return false
        end
    end

    return true
end

end
