struct GridIndex{T, Index <: AbstractVector{Int}} <: AbstractIndex
    x0::T # x value of lower edge of lowest division in x
    y0::T
    spacing::T
    n_x::Int # number of grid squares in x
    n_y::Int
    start_indices::Vector{Int} # the element of index which represents the first (lowest) point in each grid element,  with grid elements ordered from x0, y0 ... xMax,y0 ... x0,y0+1 ...etc
    index::Index # the order of the point cloud points required to match start_indices description ie pc[GridIndex.index]  will reorder the points grouping all points in the same grid element together as well as in order of z value
end

Base.summary(i::GridIndex) = "GridIndex ($(i.n_x)×$(i.n_y) cells of spacing $(i.spacing))"

Base.@propagate_inbounds Base.getindex(g::GridIndex, i, j) = g[i+g.n_x*(j-1)]
Base.@propagate_inbounds Base.getindex(g::GridIndex, i) = (g.start_indices[i]):(g.start_indices[i+1]-1)
# eg. to get all points in first grid element ie lowest x and lowest y use pc[pc.position.index.index[pc.position.index[1,1]]]

function AcceleratedArrays.accelerate(points::AbstractVector{<:StaticVector{3, T}}, ::Type{GridIndex}; spacing = one(T)) where {T}
    grid = GridIndex(points; spacing = spacing)
    return AcceleratedArray(points, grid)
end

function AcceleratedArrays.accelerate!(points::AbstractVector{<:StaticVector{3, T}}, ::Type{GridIndex}; spacing = one(T), perm = Ref{Vector{Int}}()) where {T}
    grid = GridIndex(points; spacing = spacing)
    perm[] = grid.index
    permute!(points, grid.index)
    index = keys(points)
    return AcceleratedArray(points, GridIndex(grid.x0, grid.y0, grid.spacing, grid.n_x, grid.n_y, grid.start_indices, index))
end

"""
    GridIndex(points; spacing = 1)

An acceleration index for 3D points using a 2D grid. Inside each grid cell, the points are
ordered in the third dimension. Suitable for "two-and-a-half dimensional data" such as point
clouds on the surface of the Earth, and where the query radius is close to the grid
`spacing`.

Note that `GridIndex`s are typically created by the `accelerate` and `accelerate!`
functions (which also accept the `spacing` keyword argument. In the latter case, the points
are re-ordered to increase query speed.
"""
function GridIndex(points::AbstractVector{<:StaticVector{3, T}}; spacing = one(T)) where {T}
    grid_spacing = convert(T, spacing)
    n_points = length(points)
    if n_points < 1
        return GridIndex(zero(T), zero(T), grid_spacing, 1, 1, [1, 1], Int[])
    end
    if length(points[1]) != 3
        error("Points must be of length 3")
    end

    # Initialize some grid parameters
    (xmin, xmax, ymin, ymax) = bounds(points)
    xmin = prevfloat(xmin)
    ymin = prevfloat(ymin)

    n_x = Int(cld(xmax - xmin, grid_spacing))
    n_y = Int(cld(ymax - ymin, grid_spacing))

    # Pass over the points to construct the grid
    cells = Matrix{Vector{Int}}(undef, n_x, n_y)
    for i = 1:(n_x*n_y)
        @inbounds cells[i] = Vector{Int}()
    end

    @inbounds for i ∈ 1:n_points
        p = points[i]
        cell_x = Int(cld(p[1] - xmin, grid_spacing))
        cell_y = Int(cld(p[2] - ymin, grid_spacing))
        push!(cells[cell_x, cell_y], i)
    end

    # Create the permutation and sort the cells
    permutation = Vector{Int}(undef, n_points)
    start_indices = Vector{Int}(undef, n_x::Int*n_y::Int+1)
    k = 1
    for iy ∈ 1:n_y
        for ix ∈ 1:n_x
            @inbounds start_indices[ix + n_x*(iy-1)] = k
            @inbounds cell = cells[ix,iy]
            sort!(cell, Base.Sort.QuickSortAlg(), Base.Order.By(i -> (Base.@_inline_meta; @inbounds return points[i][3])))
            k_end = k + length(cell) - 1
            @inbounds permutation[k:k_end] = cell
            k = k_end + 1
        end
    end
    @inbounds start_indices[end] = n_points+1

    grid = GridIndex(xmin, ymin, grid_spacing, n_x, n_y, start_indices, permutation)
end

function bounds(points::AbstractVector)
    PointType = eltype(points)
    T = eltype(PointType)
    xmin = typemax(T)
    xmax = typemin(T)
    ymin = typemax(T)
    ymax = typemin(T)

    @inbounds for p ∈ points
        @inbounds x = p[1]
        if x < xmin
            xmin = x
        end
        if x > xmax
            xmax = x
        end

        @inbounds y = p[2]
        if y < ymin
            ymin = y
        end
        if y > ymax
            ymax = y
        end
    end

    return (xmin, xmax, ymin, ymax)
end

function cells_bbox(grid::GridIndex, region::BoundingBox)
    # Find the rough rectangle of intersecting cells
    # Axis-aligned bounding box:
    xmin = region.xmin
    xmax = region.xmax
    ymin = region.ymin
    ymax = region.ymax

    # Make a square grid
    gx_min = Int(cld(xmin - grid.x0, grid.spacing))
    if gx_min < 1
        gx_min = 1
    end
    gx_max = Int(cld(xmax - grid.x0, grid.spacing))
    if gx_max > grid.n_x
        gx_max = grid.n_x
    end

    gy_min = Int(cld(ymin - grid.y0, grid.spacing))
    if gy_min < 1
        gy_min = 1
    end
    gy_max = Int(cld(ymax - grid.y0, grid.spacing))
    if gy_max > grid.n_y
        gy_max = grid.n_y
    end

    return (gx_min, gx_max, gy_min, gy_max)
end

# Accelerations

function Base.findall(pred::Base.Fix2{typeof(in), <:AbstractRegion}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex})
    indices = Vector{Int}()
    findall!(indices, pred, points)
    return indices
end

function findall!(indices::Vector{<:Integer}, pred::Base.Fix2{typeof(in), <:AbstractRegion}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex})
    # First find the relevant grid cells
    grid = points.index
    bbox = boundingbox(pred.x)
    (gx_min, gx_max, gy_min, gy_max) = cells_bbox(grid, bbox)

    # Now we'll iterate over these cells and fill the indices
    @inbounds for gy ∈ gy_min:gy_max
        for gx ∈ gx_min:gx_max
            for i ∈ grid[gx, gy]
                # Could use min/max height information here, as in old code
                j = grid.index[i]
                p = points[j]

                if p[3] > bbox.zmax
                    break
                end

                if pred(p)
                    push!(indices, j)
                end
            end
        end
    end

    return indices
end

function Base.count(pred::Base.Fix2{typeof(in), <:AbstractRegion}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex})
    # First find the relevant grid cells
    grid = points.index
    (gx_min, gx_max, gy_min, gy_max) = cells_bbox(grid, boundingbox(pred.x))

    out = 0

    # Now we'll iterate over these cells and fill the indices
    @inbounds for gy ∈ gy_min:gy_max
        for gx ∈ gx_min:gx_max
            for i ∈ grid[gx, gy]
                # Could use min/max height information here, as in commented code below
                j = grid.index[i]
                p = points[j]

                out += pred(p)
            end
        end
    end

    return out
end

"""
function containsMoreThanN(pred::Base.Fix2{typeof(in), <:AbstractRegion}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, N::Int64)
    Returns true if more than N points are returned as true by the predicate pred.
"""
function containsMoreThanN(pred::Base.Fix2{typeof(in), <:AbstractRegion}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, N::Int64)
    # First find the relevant grid cells
    grid = points.index
    (gx_min, gx_max, gy_min, gy_max) = cells_bbox(grid, boundingbox(pred.x))

    count = 0

    # Now we'll iterate over these cells and fill the indices
    @inbounds for gy ∈ gy_min:gy_max
        for gx ∈ gx_min:gx_max
            for i ∈ grid[gx, gy]
                j = grid.index[i]
                p = points[j]

                count += pred(p)
                if count > N
                    return true
                end
            end
        end
    end

    return false
end

function Base.filter(pred::Base.Fix2{typeof(in), <:AbstractRegion}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex})
    # First find the relevant grid cells
    grid = points.index
    (gx_min, gx_max, gy_min, gy_max) = cells_bbox(grid, boundingbox(pred.x))

    out = empty(points.parent)

    # Now we'll iterate over these cells and fill the indices
    @inbounds for gy ∈ gy_min:gy_max
        for gx ∈ gx_min:gx_max
            for i ∈ grid[gx, gy]
                # Could use min/max height information here, as in commented code below
                j = grid.index[i]
                p = points[j]

                if pred(p)
                    push!(out, p)
                end
            end
        end
    end

    return out
end

"""
findIndiciesClose2Lines(lines::AbstractVector{<:Line}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, d)

    Returns all point cloud indices inside selected grid cells, where the selected grid cells are thoes whose centre 
    are within a distance d of any of the given lines (note a 2D distance in the first 2 dimensions (ie x,y) is used).
"""
function findIndiciesClose2Lines(lines::AbstractVector{<:Line}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, d::Float64)
    grid = points.index
    lines = [convert2d(line) for line in lines]

    # Map grid cells with cell's center point
    cell_center_points = Dict()
    half_spacing = grid.spacing / 2.
    for i in 1:grid.n_x
        for j in 1:grid.n_y
            cell_center_point = SVector(grid.x0 + grid.spacing * (i - 1) + half_spacing, grid.y0 + grid.spacing * (j - 1) + half_spacing)
            cell_center_points[i, j] = cell_center_point
        end
    end

    # Find cells of center point that are in the distance d to the lines
    cell_center_point_indexs = Dict()
    for line in lines
        for i in 1:grid.n_x
            for j in 1:grid.n_y
                if !haskey(cell_center_point_indexs, (i, j))
                    cell_center_point = cell_center_points[i, j]
                    dist = distance(line, cell_center_point)
                    if dist <= d
                        cell_center_point_indexs[i, j] = nothing
                    end
                end
            end
        end
    end

    # Return indices of original point cloud of the selected cells
    inbound_points_indices = reduce(vcat, [grid.index[grid[i, j]] for (i, j) in keys(cell_center_point_indexs)]; init = Int64[])

    return inbound_points_indices
end

findIndiciesClose2Lines(points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, lines::AbstractVector{<:Line}, d::Float64) = findIndiciesClose2Lines(lines, points, d)
findIndiciesClose2Lines(points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, d::Float64, lines::AbstractVector{<:Line}) = findIndiciesClose2Lines(lines, points, d)
findIndiciesClose2Lines(lines::AbstractVector{<:Line}, d::Float64, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}) = findIndiciesClose2Lines(lines, points, d)
findIndiciesClose2Lines(d::Float64, lines::AbstractVector{<:Line}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}) = findIndiciesClose2Lines(lines, points, d)
findIndiciesClose2Lines(d::Float64, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, lines::AbstractVector{<:Line}) = findIndiciesClose2Lines(lines, points, d)

"""
Find the closest point to pos amongst all the points inside grid cell gx,gy of points. zMin and zMax provide vertical bounds within which to limit the search, best is the currenlty closest point from any potential earlier search.
minDist² is the square of the minimium distance away considerd "close enough", while maxDist² is the currenlty found closest point's distance (from previous search eg. of another cell). excLudeStart and excludeEnd allow one to only consider 
a rangge of points inside the points vector.  
"""
function closestPointInCell(pos::SVector{3,Float64}, grid::GridIndex, gx::Int64, gy::Int64, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, zMin::Float64, zMax::Float64, minDist²::Float64, maxDist²::Float64, best::Int64, excludeStart::Int64, excludeEnd::Int64)
    for i ∈ grid[gx, gy]
        j = grid.index[i]
        if excludeStart <= j <= excludeEnd
            continue
        end
        p = points[j]
        if p[3] < zMin
            continue
        end
        if p[3] > zMax
            break
        end
        diff = pos .- p
        dist² = dot(diff, diff)
        if dist² < minDist²
            return true, j, dist²
        end
        if dist² < maxDist²
            maxDist² = dist²
            best = j
        end
    end
    return false, best, maxDist²
end

"""
    findClosestPointIndex(pos::SVector{3,Float64}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, minDist::Float64, maxDist::Float64, excludeStart::Int64, excludeEnd::Int64)
  
    Looks for the closest point to pos inside the indexed array "points" stops once it finds a point within  "minDist" and will not look further than "maxDist". 
     Points between indicies excludeStart and excludeEnd are not included as part of the search.
    Returns the point index and the distance to it, will return 0 as index if nothing found. Note this code has been optimised for a particular type of point cloud
    (typical Reigl VQ 1560ii aerial capture over a combination of farmland and forest, indexed at 1m resolution) other optimisation methods could be better for point clouds of different 
    densities or index resolutions. The position pos is not assumed to exist inside points. 
"""
function findClosestPointIndex(pos::SVector{3,Float64}, points::AcceleratedArray{<:Any, <:Any, <:Any, <:GridIndex}, minDist::Float64, maxDist::Float64, excludeStart::Int64, excludeEnd::Int64)
    grid = points.index
    bbox = pad(boundingbox(pos), maxDist)
    bounds = [bbox.xmax, bbox.xmin, bbox.ymax, bbox.ymin, bbox.zmax, bbox.zmin]
    maxDist² = maxDist * maxDist
    best = 0
    minDist² = minDist * minDist

    # First look inside the index cell to which pos lies
    firstCellX = max(min(Int(cld(pos[1] - grid.x0, grid.spacing)), grid.n_x), 1)
    firstCellY = max(min(Int(cld(pos[2] - grid.y0, grid.spacing)), grid.n_y), 1)
    found, best, maxDist² =closestPointInCell(pos, grid, firstCellX, firstCellY, points, bounds[6], bounds[5], minDist², maxDist², best, excludeStart, excludeEnd)
    if found
        return best, sqrt(maxDist²)
    end

    # Find distances to edges of first cell, if all further away than best match found so far return
    dist2WEdge = (pos[1] - grid.x0)% grid.spacing
    dist2EEdge = grid.spacing - dist2WEdge
    dist2SEdge = (pos[2] - grid.y0)% grid.spacing
    dist2NEdge = grid.spacing - dist2SEdge
    dist2Edges = [dist2NEdge, dist2EEdge, dist2SEdge, dist2WEdge]
    closest = argmin(dist2Edges)
    distToEdge = dist2Edges[closest]
    if maxDist² < distToEdge * distToEdge
        return best, sqrt(maxDist²)
    end

    # check second closest cell to pos
    secondCellX = firstCellX
    secondCellY = firstCellY
    if closest == 1
        secondCellY += 1
    elseif closest == 2
        secondCellX += 1
    elseif closest == 3
        secondCellY += -1
    else
        secondCellX += -1
    end
    (gx_min, gx_max, gy_min, gy_max) = RoamesGeometry.cells_bbox(grid, bbox)
    if (gx_min < secondCellX < gx_max) && (gy_min < secondCellY < gy_max) 
        found, best, maxDist² = closestPointInCell(pos, grid, secondCellX, secondCellY, points, bounds[6], bounds[5], minDist², maxDist², best, excludeStart, excludeEnd)
        if found
            return best, sqrt(maxDist²)
        end
        dist2Edges[closest] = Inf
        closest = argmin(dist2Edges)
        distToEdge = dist2Edges[closest]
        if maxDist² < distToEdge * distToEdge
            return best, sqrt(maxDist²)
        end
    end

    # Loop through all remaining cells
    @inbounds for gy ∈ gy_min:gy_max
        for gx ∈ gx_min:gx_max

            if gx == firstCellX && gy == firstCellY || gx == secondCellX && gy == secondCellY 
                continue
            end
            
            found, best, maxDist² = closestPointInCell(pos, grid, gx, gy, points, bounds[6], bounds[5], minDist², maxDist², best, excludeStart, excludeEnd)
            if found
                return best, sqrt(maxDist²)
            end
        end
    end
    return best, sqrt(maxDist²)
end