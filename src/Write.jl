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
        println("  patch $i (start: $(patchStarts[i]) size: $(patchSizes[i])) name: $(patchNames[i])")
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
    ownerFile = open(joinpath(dir, "owner"), "w")
    neighbourFile = open(joinpath(dir, "neighbour"), "w")
    pointsFile = open(joinpath(dir, "points"), "w")
    facesFile = open(joinpath(dir, "faces"), "w")
    boundaryFile = open(joinpath(dir, "boundary"), "w")

    write_header(ownerFile, "labelList", "constant/polyMesh", "owner")
    write(ownerFile, "$(size(owner, 1))\n") 
    write(ownerFile, "(\n")
    for i in 1:size(owner, 1)
        write(ownerFile, "$(owner[i])\n")
    end
    write(ownerFile, ")\n")
    close(ownerFile)

    write_header(neighbourFile, "labelList", "constant/polyMesh", "neighbour")
    write(neighbourFile, "$(size(neighbour, 1))\n") 
    write(neighbourFile, "(\n")
    for i in 1:size(neighbour, 1)
        write(neighbourFile, "$(neighbour[i])\n")
    end
    write(neighbourFile, ")\n")
    close(neighbourFile)

    write_header(pointsFile, "vectorField", "constant/polyMesh", "points")
    write(pointsFile, "$(size(points, 1))\n") 
    write(pointsFile, "(\n")
    for i in 1:size(points, 1)
        write(pointsFile, "($(points[i][1]) $(points[i][2]) $(points[i][3]))\n")
    end
    write(pointsFile, ")\n")
    close(pointsFile)

    write_header(facesFile, "faceList", "constant/polyMesh", "faces")
    write(facesFile, "$(size(faces, 1))\n") 
    write(facesFile, "(\n")
    for i in 1:size(faces, 1)
        write(facesFile, "4($(faces[i][1]) $(faces[i][2]) $(faces[i][3]) $(faces[i][4]))\n")
    end
    write(facesFile, ")\n")
   close(facesFile)

    nPatches = length(patchDicts)
    write_header(boundaryFile, "polyBoundaryMesh", "constant/polyMesh", "boundary")
    write(boundaryFile, "$(nPatches)\n") 
    write(boundaryFile, "(\n")
    for i in 1:nPatches
        write(boundaryFile, "    $(patchDicts[i]["name"])\n")
        write(boundaryFile, "    {\n")
        write(boundaryFile, "    type      $(patchDicts[i]["type"]);\n")
        for (key, val) in patchDicts[i]
            if !(key in ["name", "type", "faces"])
                write(boundaryFile, "    $(key) $(val);\n")
            end
        end
        write(boundaryFile, "    nFaces    $(patchSizes[i]);\n")
        write(boundaryFile, "    startFace $(patchStarts[i]);\n")
        write(boundaryFile, "    }\n")

    end

    write(boundaryFile, ")\n")
    close(boundaryFile)
end
