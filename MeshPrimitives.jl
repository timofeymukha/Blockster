module MeshPrimitives

using FixedSizeArrays

export point, cell


" Type for a point in 3d space. An alias for a vector of doubles of size 3"
typealias point Vec{3, Float64}

typealias cell Vec{8, Int64 }

end
