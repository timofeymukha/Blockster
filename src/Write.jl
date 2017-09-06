export write_header, write_mesh_information, write_mesh

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

function write_mesh_information(
    nPoints,
    nCells,
    nFaces,
    nInternalFaces,
    patchNames,
    patchStarts,
    patchSizes
)
    println("----------------")
    println("Mesh Information")
    println("----------------")
    println("  nPoints $(nPoints)")
    println("  nCells $(nCells)")
    println("  nFaces $(nFaces)")
    println("  nInternalFaces $(nInternalFaces)")

    println("----------------")
    println("Patches")
    println("----------------")

    for i in 1:length(patchSizes)
        println("  patch $(i-1) (start: $(patchStarts[i]) size: $(patchSizes[i])) name: $(patchNames[i])")
    end

end

function write_mesh(
    dir,
    owner,
    neighbour,
    points,
    faces,
    patchStarts,
    patchSizes,
    patchDicts
)

    @inbounds begin
    ownerFile = open(joinpath(dir, "owner"), "w")
    neighbourFile = open(joinpath(dir, "neighbour"), "w")
    pointsFile = open(joinpath(dir, "points"), "w")
    facesFile = open(joinpath(dir, "faces"), "w")
    boundaryFile = open(joinpath(dir, "boundary"), "w")

    write_header(ownerFile, "labelList", "constant/polyMesh", "owner")
    println(ownerFile, size(owner, 1)) 
    println(ownerFile, "(")
    for i in 1:size(owner, 1)
        println(ownerFile, owner[i])
    end
    println(ownerFile, ")")
    close(ownerFile)

    write_header(neighbourFile, "labelList", "constant/polyMesh", "neighbour")
    println(neighbourFile, length(neighbour)) 
    write(neighbourFile, "(\n")
    for i in 1:size(neighbour, 1)
        println(neighbourFile, neighbour[i])
    end
    println(neighbourFile, ")")
    close(neighbourFile)

    write_header(pointsFile, "vectorField", "constant/polyMesh", "points")
    println(pointsFile, length(points)) 
    println(pointsFile, "(")
    for i in 1:size(points, 1)
        println(pointsFile, "(", points[i][1], " ", points[i][2], " ", points[i][3], ")")
    end
    println(pointsFile, ")")
    close(pointsFile)

    write_header(facesFile, "faceList", "constant/polyMesh", "faces")
    println(facesFile, length(faces)) 
    println(facesFile, "(")
    for i in 1:size(faces, 1)
        println(facesFile, "4(", faces[i][1], " ", faces[i][2], " ", faces[i][3], " ", faces[i][4], ")")
    end
    println(facesFile, ")")
    close(facesFile)

    nPatches = length(patchDicts)
    write_header(boundaryFile, "polyBoundaryMesh", "constant/polyMesh", "boundary")
    println(boundaryFile, nPatches) 
    println(boundaryFile, "(")
    for i in 1:nPatches
        println(boundaryFile, "    $(patchDicts[i]["name"])")
        println(boundaryFile, "    {")
        println(boundaryFile, "    type      $(patchDicts[i]["type"]);")
        for (key, val) in patchDicts[i]
            if !(key in ["name", "type", "faces"])
                println(boundaryFile, "    $(key) $(val);")
            end
        end
        write(boundaryFile, "    nFaces    $(patchSizes[i]);\n")
        write(boundaryFile, "    startFace $(patchStarts[i]);\n")
        println(boundaryFile, "    }")

    end

    println(boundaryFile, ")")
    close(boundaryFile)
    end #inbounds
end
