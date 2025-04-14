struct Polygon{N, T <: Real, L <: LineString{N, T}, V <: AbstractVector{L}} <: AbstractRegion{T}
    exterior::L
    interiors::V
end

if VERSION < v"0.7"
    function Polygon(exterior::L, interiors::AbstractVector{L}) where {N, T, L <: LineString{N, T}}
        Polygon{N, T, L, typeof(interiors)}(exterior, interiors)
    end
end
function Polygon{N}(exterior::L, interiors::AbstractVector{L}) where {N, T, L <: LineString{N, T}}
    Polygon{N, T, L, typeof(interiors)}(exterior, interiors)
end
function Polygon{N,T}(exterior::L, interiors::AbstractVector{L}) where {N, T, L <: LineString{N, T}}
    Polygon{N, T, L, typeof(interiors)}(exterior, interiors)
end

Polygon(points::StaticVector{N,T}...) where {N, T<:Real} = Polygon(collect(points))
Polygon(points::AbstractVector{<:StaticVector{N,T}}) where {N,T<:Real} = Polygon{N}(LineString{N}(points))
Polygon{N}(points::AbstractVector{T}...) where {N, T<:Real} = Polygon{N}(collect(points))
Polygon{N}(points::AbstractVector{<:AbstractVector{T}}) where {N, T<:Real} = Polygon{N,T}(LineString{N,T}(points))
Polygon{N,T}(points::AbstractVector{<:Real}...) where {N, T<:Real} = Polygon{N,T}(collect(points))
Polygon{N,T}(points::AbstractVector{<:AbstractVector{<:Real}}) where {N, T<:Real} = Polygon{N,T}(LineString{N,T}(points))

Polygon(ls::LineString) = Polygon(ls, Vector{typeof(ls)}())
Polygon{N}(ls::LineString{N}) where {N} = Polygon{N}(ls, Vector{typeof(ls)}())
Polygon{N,T}(ls::LineString{N,T}) where {N,T<:Real} = Polygon{N,T}(ls, Vector{typeof(ls)}())

convert2d(p::Polygon{2}) = p
convert2d(p::Polygon{3}) = Polygon(convert2d(p.exterior), convert2d.(p.interiors))
convert3d(p::Polygon{2}) = Polygon(convert3d(p.exterior), convert3d.(p.interiors))
convert3d(p::Polygon{2}, z::Real) = Polygon(convert3d(p.exterior, z), convert3d.(p.interiors, z))
convert3d(p::Polygon{3}) = p
convert3d(p::Polygon{3}, z::Real) = Polygon(convert3d(p.exterior, z), convert3d.(p.interiors, z))

# Note: This doesn't compare the ordering of points
function Base.:(==)(p1::Polygon{N}, p2::Polygon{N}) where N
    p1.exterior == p2.exterior && p1.interiors == p2.interiors
end

function Base.isequal(p1::Polygon{N}, p2::Polygon{N}) where N
    isequal(p1.exterior, p2.exterior) && isequal(p1.interiors, p2.interiors)
end

function Base.hash(p::Polygon, h::UInt)
    hash(p.exterior, hash(p.interiors, hash(UInt === UInt64 ? 0xde95b490c51c55a5 : 0x0f7345a4, h)))
end

function show(io::IO, polygon::Polygon{N}) where N
    print(io, "Polygon{$N}([")
    for i in 1:length(polygon.exterior.points)
        print(io, polygon.exterior.points[i])
        if i < length(polygon.exterior.points)
            print(io, ", ")
        end
    end
    print(io, "]")
    for ls in polygon.interiors
    	print(io, ", [")
    	for i in 1:length(ls.points)
	        print(io, ls.points[i])
	        if i < length(ls.points)
	            print(io, ", ")
	        end
	    end
    	print(io, "]")
    end
    print(io, ")")
end

function lines(p::Polygon)
    out = lines(p.exterior)
    for ls in p.interiors
        append!(out, lines(ls))
    end
    return out
end

function area(p::Polygon{2, T}) where {T}
    return area(p.exterior) + mapreduce(area, +, p.interiors; init = zero(T))
end

function in(p::StaticVector{2, <:Real}, multiPolygon::Vector{P} where P<:Polygon)
    return any(poly->p in poly, multiPolygon)
end

function in(p::StaticVector{2, <:Real}, polygon::Polygon{2})
    # Winding number algorith, for example read
    # http://geomalgorithms.com/a03-_inclusion.html

    # Count number of times an edge of the polygon crosses a ray from
    # p to infinity (doesn't matter which direction but we choose +x direction)
    wn = winding_number(p, polygon.exterior)
    for linestring in polygon.interiors
        wn += winding_number(p, linestring)
    end

    return wn != 0 # Valid polygons can be oriented clockwise or anticlockwise
end

intersects(p::StaticVector{2, <:Real}, polygon::Polygon{2}) = p ∈ polygon
intersects(polygon::Polygon{2}, p::StaticVector{2, <:Real}) = p ∈ polygon

function intersects(l::Line{2}, polygon::Polygon{2})
    # The polygon is a surface. Either `l` intersects the edge of
    # the polygon or it is entirely *inside* the polygon (or else it
    # doesn't intersect at all)
    if intersects(l, polygon.exterior)
        return true
    end

    for linestring in polygon.interiors
        if intersects(l, linestring)
            return true
        end
    end

    # Either both endpoints are inside, or outside, so just test one
    return l.p1 ∈ polygon
end
intersects(polygon::Polygon{2}, l::Line{2}) = intersects(l, polygon)

function intersects(ls::LineString{2}, polygon::Polygon{2})
    # The polygon is a surface. Either `ls` intersects the edge of
    # the polygon, or it is entirely *inside* the polygon (or else it
    # doesn't intersect at all)
    if intersects(ls, polygon.exterior)
        return true
    end

    for linestring in polygon.interiors
        if intersects(ls, linestring)
            return true
        end
    end

    for p in ls.points
        if p ∉ polygon
            return false
        end
    end

    # `ls` must be entirely inside `polygon`
    return true
end
intersects(polygon::Polygon{2}, ls::LineString{2}) = intersects(ls, polygon)

function intersects(p1::Polygon{2}, p2::Polygon{2})
    # If any line (interior OR exterior) intersects with any other line,
    # then the two polygons intersect. If not, check if one polygon
    # is entirely inside the other
    
    # Check if edges intersect
    if intersects(p1.exterior, p2)
        return true
    end

    for linestring in p1.interiors
        if intersects(linestring, p2)
            return true
        end
    end

    # One polygon might be entirely inside another
    # Since there can't be multple "onion" rings in any
    # polygon (it can have holes but is otherwise contiguous)
    # we only need to check the exteriors
    if all(p -> p in p2, p1.exterior.points)
        return true
    end
    if all(p -> p in p1, p2.exterior.points)
        return true
    end

    # Or else they do not intersect
    return false
end

function intersects(p1::Polygon{2}, p2::Vector{Polygon})
    return any(pp-> intersects(pp, p1) ,p2)
end
intersects(p2::Vector{Polygon}, p1::Polygon{2}) = intersects(p1, p2)

function issimple(polygon::Polygon)
    return issimple(polygon.exterior) && all(p->issimple(p), polygon.interiors)
end

function issimple(multiPolygon::Vector{P} where P <: Polygon)
    return all(p->issimple(p), multiPolygon)
end

"""
function RoamesPolyToLibGEOSPoly(RoamesPolygon::Polygon)

Convert a Roames Polygon to a LibGEOS Polygon
Note a 3D polygon will become 2D after conversion.
"""
function RoamesPolyToLibGEOSPoly(RoamesPolygon::Polygon)
    exteriorRingCoords = map(p->p[1:2], RoamesPolygon.exterior.points)

    # Create the exterior ring in LibGEOS format
    exteriorRing = LibGEOS.createLinearRing(exteriorRingCoords)
    if exteriorRing == C_NULL
        error("Failed to create exterior ring")
    end

    # Handle interior rings (holes) if they exist
    numInteriorRings = length(RoamesPolygon.interiors)
    interiorRings = Vector{Ptr{LibGEOS.GEOSGeom}}(undef, numInteriorRings)

    for i in 1:numInteriorRings
        interiorRingCoords = map(p->[p[1:2]], RoamesPolygon.interiors[i].points)
        interiorRings[i] = LibGEOS.createLinearRing(interiorRingCoords)

        if interiorRings[i] == C_NULL
            error("Failed to create Lib Geos interior ring $i from Roames Polygon $(RoamesPolygon) ")
        end
    end

    # Create the LibGEOS polygon
    LibGEOSPolygon = LibGEOS.createPolygon(exteriorRing,  interiorRings)

    if LibGEOSPolygon == C_NULL
        error("Failed to create polygon")
    end

    geom = LibGEOS.geomFromGEOS(LibGEOSPolygon)
    if !LibGEOS.isValid(geom)
        error(LibGEOS.isValidReason(geom))
    end

    return geom
end

"""
function RoamesPolyToLibGEOSPoly(RoamesPolygon::Vector{Polygon})

Convert a vector of Roames Polygons to a LibGEOS MultiPolygon
Note a 3D polygon will become 2D after conversion.
"""
function RoamesPolyToLibGEOSPoly(RoamesPolygon::Vector{P} where P < Polygon)
    return LibGEOS.MultiPolygon(RoamesPolyToLibGEOSPoly.(RoamesPolygon))
end

"""
function LibGEOSPolyToRoamesPoly(LibGEOSPolygon::LibGEOS.Polygon)

Convert a libgeos Polygon to a Roames Polygon
"""
function LibGEOSPolyToRoamesPoly(LibGEOSPolygon::LibGEOS.Polygon)
    # Get the exterior ring
    exteriorRingGeos = LibGEOS.exteriorRing(LibGEOSPolygon)
    if exteriorRingGeos == C_NULL
        error("Polygon has no exterior ring")
    end

    # Get the coordinate sequence of the exterior ring
    coordSeq = LibGEOS.getCoordSeq(exteriorRingGeos)
    if coordSeq == C_NULL
        error("Failed to get coordinate sequence for exterior ring")
    end

    exteriorCoords = LibGEOS.getCoordinates(coordSeq)
    exteriorRingCoords = [ SVector(x, y) for (x, y) in exteriorCoords]
    outer = LineString(exteriorRingCoords)

    # inter ring goes here
    numInteriorRings = LibGEOS.numInteriorRings(LibGEOSPolygon)
    inners = LineString{2, Float64, Vector{RoamesGeometry.SVector{2, Float64}}}[]
    for i in 1:numInteriorRings
        interiorRingGeos = LibGEOS.interiorRing(LibGEOSPolygon,i)
        coordSeq = LibGEOS.getCoordSeq(interiorRingGeos)
        interiorCoords = LibGEOS.getCoordinates(coordSeq)
        interiorRingCoords = [ SVector(x, y) for (x, y) in interiorCoords]
        push!(inners,LineString(interiorRingCoords))
    end

    Rpolygon = Polygon(outer, inners)
    if Rpolygon == nothing
        error("Failed to convert to RoamesGeometry.Polygon")
    end
    return Rpolygon
end

"""
function LibGEOSPolyToRoamesPoly(LibGEOSPolygon::LibGEOS.MultiPolygon)

Convert a libgeos MultiPolygon to a RoamesGeometry Vector of Polygons
"""
function LibGEOSPolyToRoamesPoly(LibGEOSPolygon::LibGEOS.MultiPolygon)
    return map(p->LibGEOSPolyToRoamesPoly(RoamesGeometry.LibGEOS.getgeom(LibGEOSPolygon, p)), 1:LibGEOS.ngeom(LibGEOSPolygon))
end
