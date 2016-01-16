module Blocks

using FixedSizeArrays
using Edges
using MeshPrimitives

export Block, make_block_edges!, setedge!, create_points!, create_cells!,
       point_index, npoints, ncells 

" Type that defines a single block of the multi-block mesh."
type Block
    vertexLabels::Vector{Int64}
    vertices::Vector{point}
    points::Vector{point}
    cells::Array{Int64,2}
    edgePoints::Vector{Vector{point}}
    edgeWeights::Vector{Vector{Float64}}
    curvedEdges::Vector{CurvedEdge}
    nCells::Vec{3, Int64}
end

Block() = Block(zeros(8),
                Vector{point}(8),
                Vector{point}(0),
                zeros(1,8),
                Vector{Vector{point}}(12),
                Vector{Vector{Float64}}(12),
                Vector{CurvedEdge}(0),
                zeros(3))



function make_block_edges!(block::Block)
    nX = block.nCells[1];
    nY = block.nCells[2];
    nZ = block.nCells[3];

    # These edges correspond to the "hex" cellModel

    # X-direction
    setedge!(block, 1,  1, 2, nX);
    setedge!(block, 2,  4, 3, nX);
    setedge!(block, 3,  8, 7, nX);
    setedge!(block, 4,  5, 6, nX);

    # Y-direction
    setedge!(block, 5,  1, 4, nY);
    setedge!(block, 6,  2, 3, nY);
    setedge!(block, 7,  6, 7, nY);
    setedge!(block, 8,  5, 8, nY);

    # Z-direction
    setedge!(block, 9,  1, 5, nZ);
    setedge!(block, 10, 2, 6, nZ);
    setedge!(block, 11, 3, 7, nZ);
    setedge!(block, 12, 4, 8, nZ);
end

function setedge!(block::Block, edgeI::Int64, startVertex::Int64,
                 endVertex::Int64, nDivisions::Int64)
    # Check that start and end vertices are different
    if startVertex == endVertex
        warn("setedge!() : start and end vertices are the same.")
    end

    # Loop through the curved edges and compare to given edge
    for i in 1:size(block.curvedEdges, 1)

        curvedEdge = block.curvedEdges[i]

        cmp = compare(curvedEdge,
                    block.vertices[startVertex],
                    block.vertices[endVertex])

        if cmp > 0
            # Divide the line
        else
            # Divide the line

            # Found and divided the edge, exit
            return
        end
    end
    # The edge is a straight line

    edge = LineEdge(block.vertices, startVertex, endVertex)

    (block.edgePoints[edgeI],
     block.edgeWeights[edgeI]) = line_divide(edge, nDivisions)
end

function point_index(block::Block, i::Int64, j::Int64, k::Int64)
    # Check bounds
    if i < 1 || j < 1 || k < 1
        error("point_index(): index less then 1")
    elseif i > block.nCells[1] + 1 ||
           j > block.nCells[2] + 1 ||
           k > block.nCells[3] + 1
        error("point_index(): index out of bounds")
    end
        
    return i + (j-1)*(block.nCells[1] + 1) + 
               (k-1)*(block.nCells[1] + 1)*(block.nCells[2] + 1)
end

function npoints(block::Block)
    return (block.nCells[1] + 1)*(block.nCells[2] + 1)*(block.nCells[3] + 1)
end

function ncells(block::Block)
    return block.nCells[1]*block.nCells[2]*block.nCells[3]
end

function create_points!(block::Block)

    v000 = block.vertices[1];
    v100 = block.vertices[2];
    v110 = block.vertices[3];
    v010 = block.vertices[4];

    v001 = block.vertices[5];
    v101 = block.vertices[6];
    v111 = block.vertices[7];
    v011 = block.vertices[8];

    edgeP = block.edgePoints
    edgeW = block.edgeWeights

    block.points = resize!(block.points, npoints(block))
    
    for k in 1:block.nCells[3]+1
        for j in 1:block.nCells[2]+1
            for i in 1:block.nCells[1]+1
                pointIndex = point_index(block, i, j, k);

                edgeX1 = v000 + (v100 - v000)*edgeW[1][i];
                edgeX2 = v010 + (v110 - v010)*edgeW[2][i];
                edgeX3 = v011 + (v111 - v011)*edgeW[3][i];
                edgeX4 = v001 + (v101 - v001)*edgeW[4][i];

                edgeY1 = v000 + (v010 - v000)*edgeW[5][j];
                edgeY2 = v100 + (v110 - v100)*edgeW[6][j];
                edgeY3 = v101 + (v111 - v101)*edgeW[7][j];
                edgeY4 = v001 + (v011 - v001)*edgeW[8][j];

                edgeZ1 = v000 + (v001 - v000)*edgeW[9][k];
                edgeZ2 = v100 + (v101 - v100)*edgeW[10][k];
                edgeZ3 = v110 + (v111 - v110)*edgeW[11][k];
                edgeZ4 = v010 + (v011 - v010)*edgeW[12][k];
                
                # Calculate the importance factors for all edges

                # x-direction
                impX1 = (
                        (1 - edgeW[1][i])*(1 - edgeW[5][j])*(1 - edgeW[9][k])
                        + edgeW[1][i]*(1 - edgeW[6][j])*(1 - edgeW[10][k])
                        )
                impX2 = (
                        (1 - edgeW[2][i])*edgeW[5][j]*(1 - edgeW[12][k])
                        + edgeW[2][i]*edgeW[6][j]*(1 - edgeW[11][k])
                        )
                impX3 = (
                        (1 - edgeW[3][i])*edgeW[8][j]*edgeW[12][k]
                        + edgeW[3][i]*edgeW[7][j]*edgeW[11][k]
                        )
                impX4 = (
                        (1 - edgeW[4][i])*(1 - edgeW[8][j])*edgeW[9][k]
                        + edgeW[4][i]*(1 - edgeW[7][j])*edgeW[10][k]
                        )
                magImpX = impX1 + impX2 + impX3 + impX4
                impX1 /= magImpX
                impX2 /= magImpX
                impX3 /= magImpX
                impX4 /= magImpX


                # y-direction
                impY1 = (
                        (1 - edgeW[5][j])*(1 - edgeW[1][i])*(1 - edgeW[9][k])
                        + edgeW[5][j]*(1 - edgeW[2][i])*(1 - edgeW[12][k])
                        )

                impY2 = (
                        (1 - edgeW[6][j])*edgeW[1][i]*(1 - edgeW[10][k])
                        + edgeW[6][j]*edgeW[2][i]*(1 - edgeW[11][k])
                        )

                impY3 = (
                        (1 - edgeW[7][j])*edgeW[4][i]*edgeW[10][k]
                        + edgeW[7][j]*edgeW[3][i]*edgeW[11][k]
                        )

                impY4 = (
                        (1 - edgeW[8][j])*(1 - edgeW[4][i])*edgeW[9][k]
                        + edgeW[8][j]*(1 - edgeW[3][i])*edgeW[12][k]
                        )

                magImpY = impY1 + impY2 + impY3 + impY4

                impY1 /= magImpY
                impY2 /= magImpY
                impY3 /= magImpY
                impY4 /= magImpY



                # z-direction
                impZ1 = (
                        (1 - edgeW[9][k])*(1 - edgeW[1][i])*(1 - edgeW[5][j])
                        + edgeW[9][k]*(1 - edgeW[4][i])*(1 - edgeW[8][j])
                        )

                impZ2 = (
                        (1 - edgeW[10][k])*edgeW[1][i]*(1 - edgeW[6][j])
                        + edgeW[10][k]*edgeW[4][i]*(1 - edgeW[7][j])
                        )

                impZ3 = (
                        (1 - edgeW[11][k])*edgeW[2][i]*edgeW[6][j]
                        + edgeW[11][k]*edgeW[3][i]*edgeW[7][j]
                        )

                impZ4 = (
                        (1 - edgeW[12][k])*(1 - edgeW[2][i])*edgeW[5][j]
                        + edgeW[12][k]*(1 - edgeW[3][i])*edgeW[8][j]
                        )

                magImpZ = impZ1 + impZ2 + impZ3 + impZ4;
                impZ1 /= magImpZ
                impZ2 /= magImpZ
                impZ3 /= magImpZ
                impZ4 /= magImpZ


                # Calculate the correction vectors
                corX1 = impX1*(edgeP[1][i] - edgeX1);
                corX2 = impX2*(edgeP[2][i] - edgeX2);
                corX3 = impX3*(edgeP[3][i] - edgeX3);
                corX4 = impX4*(edgeP[4][i] - edgeX4);

                corY1 = impY1*(edgeP[5][j] - edgeY1);
                corY2 = impY2*(edgeP[6][j] - edgeY2);
                corY3 = impY3*(edgeP[7][j] - edgeY3);
                corY4 = impY4*(edgeP[8][j] - edgeY4);

                corZ1 = impZ1*(edgeP[9][k] - edgeZ1);
                corZ2 = impZ2*(edgeP[10][k] - edgeZ2);
                corZ3 = impZ3*(edgeP[11][k] - edgeZ3);
                corZ4 = impZ4*(edgeP[12][k] - edgeZ4);


                # Multiply by the importance factor

                # x-direction
                edgeX1 *= impX1;
                edgeX2 *= impX2;
                edgeX3 *= impX3;
                edgeX4 *= impX4;

                # y-direction
                edgeY1 *= impY1;
                edgeY2 *= impY2;
                edgeY3 *= impY3;
                edgeY4 *= impY4;

                # z-direction
                edgeZ1 *= impZ1;
                edgeZ2 *= impZ2;
                edgeZ3 *= impZ3;
                edgeZ4 *= impZ4;


                # Add the contributions
                block.points[pointIndex] = (edgeX1 + edgeX2 + edgeX3 + edgeX4
                                           + edgeY1 + edgeY2 + edgeY3 + edgeY4
                                           + edgeZ1 + edgeZ2 + edgeZ3 + edgeZ4
                                           )/3

                block.points[pointIndex] += (corX1 + corX2 + corX3 + corX4
                                            + corY1 + corY2 + corY3 + corY4
                                            + corZ1 + corZ2 + corZ3 + corZ4)
            end
        end
    end
end

function create_cells!(block::Block)

    block.cells = zeros(ncells(block), 8)

    cellNo::Int64 = 1

    for k in 1:block.nCells[3]
        for j in 1:block.nCells[2]
            for i in 1:block.nCells[1]
                block.cells[cellNo, 1] =  point_index(block, i, j, k)
                block.cells[cellNo, 2] =  point_index(block, i+1, j, k)
                block.cells[cellNo, 3] =  point_index(block, i+1, j+1, k)
                block.cells[cellNo, 4] =  point_index(block, i, j+1, k)
                block.cells[cellNo, 5] =  point_index(block, i, j, k+1)
                block.cells[cellNo, 6] =  point_index(block, i+1, j, k+1)
                block.cells[cellNo, 7] =  point_index(block, i+1, j+1, k+1)
                block.cells[cellNo, 8] =  point_index(block, i, j+1, k+1)
            
                cellNo += 1  
            end
        end
    end
end

end
