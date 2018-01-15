replacecoordinate!(x, i::Integer, val) = (x[i] = val; x)

replacecoordinate!(x::SVector{N,T}, i::Integer, val) where {N,T} =
    SVector{N,T}(_rpc(Tuple(x), i-1, T(val)))
@inline _rpc(t, i, val) = (ifelse(i == 0, val, t[1]), _rpc(tail(t), i-1, val)...)
_rps(::Tuple{}, i, val) = ()

ipcopy!(dest, src) = copy!(dest, src)
ipcopy!(dest::SVector, src) = src

"""
    ileft, iright, iclosest = splitindexes(box)

Among the evaluated split points in `box.parent`, which correspond to
the closest points to `box`. `ileft` is the index of the point to the "left,"
`iright` to the right, and `iclosest` is whichever of these actually lies on
one of the boundaries of `box`.
"""
function splitindexes(box)
    isroot(box) && error("root box doesn't have a parent")
    p = box.parent
    cindex = box.parent_cindex
    overhang_left = p.children[1].edge < p.xvalues[1]
    offset = (1-2*overhang_left)/4 + 1
    xn = (cindex - 1)/2 + offset
    floor(Int, xn), ceil(Int, xn), round(Int, xn)
end

"""
   xmin, fmin = qfit_min((xm, fm), (x0, f0), (xp, fp))

Given three points `xm < x0 < xp` and three corresponding values `fm`, `f0`, and `fp`,
fit a quadratic. If the quadratic is convex, return the position of the minimum and
the corresponding value of the quadratic.

If the quadratic is not convex, returns a point outside the interval by an amount equal to
the separation of the nearest pair (e.g., if f is increasing, returns a value `xm - (x0-xm)`).
"""
function qfit_min(xfm, xf0, xfp)
    xm, fm = xfm
    x0, f0 = xf0
    xp, fp = xfp
    @assert(xp > x0 && x0 > xm && isfinite(xm) && isfinite(xp))
    cm = fm/((xm-x0)*(xm-xp))  # coefficients of Lagrange polynomial
    c0 = f0/((x0-xm)*(x0-xp))
    cp = fp/((xp-xm)*(xp-x0))
    csum = cm+c0+cp
    qvalue(x) = cm*(x-x0)*(x-xp) + c0*(x-xm)*(x-xp) + cp*(x-xm)*(x-x0)
    if csum > 0
        # It's convex
        xmin = (cm*(x0+xp) + c0*(xm+xp) + cp*(xm+x0))/(2*csum)
        return xmin, qvalue(xmin)
    end
    if fm == f0 == fp
        return x0, f0
    end
    @assert(f0 >= fm || f0 >= fp)
    # This should only happen at a box edge. Just return something
    # to indicate the direction of the minimum
    if f0 > fm && fp > fm
        xmin = xm - (x0 - xm)
    else
        xmin = xp + (xp - x0)
    end
    return xmin, qvalue(xmin)
end

function Base.minimum(box::Box)
    isroot(box) && error("cannot compute minimum for the root box")
    ileft, iright, iclosest = splitindexes(box)
    p = box.parent
    (ileft == 0 || iright > length(p.xvalues)) && return p.fvalues[iclosest]
    # Compute value on "internal" edge (between xvalues) by linear interpolation
    frac = abs(box.edge - p.xvalues[iclosest])/(p.xvalues[iright] - p.xvalues[ileft])
    iother = iclosest == ileft ? iright : ileft
    fc, fo = p.fvalues[iclosest], p.fvalues[iother]
    return min(fc, fc + (fo-fc)*frac)
end

## Tree traversal
function get_root(box::Box)
    while !isroot(box)
        box = box.parent
    end
    box
end

struct DepthFirstIterator{T}
    root::Box{T}
end
struct DepthFirstState{T}
    cparent::Box{T}
    cindex::Int
end

function visit_leaves(root::Box)
    DepthFirstIterator(root)
end

function Base.start(iter::DepthFirstIterator)
    find_next_leaf(iter, DepthFirstState(iter.root, 0))
end
function Base.done(iter::DepthFirstIterator, state::DepthFirstState)
    !issplit(iter.root) && return true
    iter.root == state.cparent && state.cindex > length(iter.root.children)
end
function Base.next(iter::DepthFirstIterator, state::DepthFirstState)
    box = state.cparent.children[state.cindex]
    @assert(!issplit(box))
    return (box, find_next_leaf(iter, state))
end
function find_next_leaf(iter::DepthFirstIterator, state::DepthFirstState)
    !issplit(iter.root) && return state
    box, i = state.cparent, state.cindex+1
    if i > length(box.children)
        box, i = up(box, iter.root)
    end
    if i <= length(box.children) && issplit(box.children[i])
        return find_next_leaf(iter, DepthFirstState(box.children[i], 0))
    end
    DepthFirstState(box, i)
end
function up(box, root)
    local i
    while true
        box, i = box.parent, box.parent_cindex+1
        box == root && return (box, i)
        i <= length(box.children) && break
    end
    return (box, i)
end
