
struct LineEdge <: CurvedEdge
    points::Vector{Point}
    startVertex::Int64
    endVertex::Int64
end

function position(edge::LineEdge, lambda::Float64)
    if (lambda < 0) || (lambda > 1)
        error("position() : Parameter out of range, lambda should be between 0 and 1")
    end

    return edge.points[edge.startVertex] + 
           lambda*(edge.points[edge.endVertex] - edge.points[edge.startVertex])
end

function mag(edge::LineEdge)
    return abs(edge.points[endVertex] - edge.points[startVertex]);
end

