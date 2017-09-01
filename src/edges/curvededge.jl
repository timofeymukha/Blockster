abstract  type CurvedEdge
end

function compare(
    edge::CurvedEdge,
    startVertex::Int64,
    endVertex::Int64
)
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

function compare(
    edge1::CurvedEdge,
    edge2::CurvedEdge
)
    return compare(edge1, edge2.startVertex, edge2.endVertex)
end

function line_divide(edge::CurvedEdge, nDivisions::Int, grading)

    points = Vector{Point}(nDivisions+1)
    divisions = Vector{Float64}(nDivisions+1)
    
    divisions[1] = 0
    divisions[end] = 1

    sectionStart::Float64 = divisions[1]
    sectionStartInd::Int = 2

    if (nDivisions >= length(grading))

        sectionDivs = Vector{Int}(length(grading))
        sumSectionDivs = 0
        sectionMaxDivs = 1

        for sectionI in 1:length(grading)
            
            nDivFraction = grading[sectionI][2]
            sectionDivs[sectionI] = round(nDivFraction*nDivisions + 0.5)
            sumSectionDivs += sectionDivs[sectionI]

            # Find the section with the largest number of divisions

            if nDivFraction > grading[sectionMaxDivs][2]
                sectionMaxDivs = sectionI
            end
        end


        if sumSectionDivs != nDivisions
            sectionDivs[sectionMaxDivs] += (nDivisions - sumSectionDivs)
        end

        for sectionI in 1:length(grading) 

            blockFraction = grading[sectionI][1]
            expRatio = grading[sectionI][3]
            
            sectionDiv = sectionDivs[sectionI]
            sectionEndInd::Int = sectionStartInd + sectionDiv

            if expRatio == 1
                for i = sectionStartInd:sectionEndInd-1
                    divisions[i] =
                        sectionStart +
                        blockFraction*(i - sectionStartInd + 1)/sectionDiv
                end

            else
                expFact = cell_to_cell_expansion_ratio(expRatio, sectionDiv)
                for i = sectionStartInd:sectionEndInd - 1
                    divisions[i] =
                        sectionStart +
                        blockFraction*(1 - expFact^(i - sectionStartInd + 1))/
                        (1 - expFact^sectionDiv)
                end
            end

            sectionStart = divisions[sectionEndInd - 1]
            sectionStartInd = sectionEndInd
        end
    else

        for i in 1:nDivisions
            divisions[i] = (i-1)/nDivisions 
        end
    end
    
    # Calculate te points
    for i in 1:nDivisions+1
        points[i] = position(edge, divisions[i])
    end

    return (points, divisions)
end

function cell_to_cell_expansion_ratio(expRatio::Float64, nDivisions::Int)
    return nDivisions > 1 ? expRatio^(1.0/(nDivisions - 1)) : 0.0
end


