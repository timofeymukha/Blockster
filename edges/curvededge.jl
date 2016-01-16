abstract CurvedEdge

function compare(edge::CurvedEdge,
                 startVertex::Int64,
                 endVertex::Int64)

    if edge.startVertex == startVertex &&
       edge.endVertex == endVertex
        return 1
    elseif edge.startVertex == endVertex &&
           edge.endVertex == startVertex
        return -1
    else
        return 0
    end
end

function compare(edge1::CurvedEdge, edge2::CurvedEdge)
    return compare(edge1, edge2.startVertex, edge2.endVertex)
end

function line_divide(edge::CurvedEdge, nDivisions::Int64)
    points = Vector{Point}(nDivisions+1)
    divisions = Vector{Float64}(nDivisions+1)
    
    divisions[end] = 1
    # Should treat general edge division with grading.

    for i in 1:nDivisions
        divisions[i] = (i-1)/nDivisions 
    end
    
    for i in 1:nDivisions+1
        points[i] = position(edge, divisions[i])
    end

    return (points, divisions)
end
