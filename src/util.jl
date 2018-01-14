replacecoordinate!(x, i::Integer, val) = (x[i] = val; x)

replacecoordinate!(x::SVector{N,T}, i::Integer, val) where {N,T} =
    SVector{N,T}(_rpc(Tuple(x), i-1, T(val)))
@inline _rpc(t, i, val) = (ifelse(i == 0, val, t[1]), _rpc(tail(t), i-1, val)...)
_rps(::Tuple{}, i, val) = ()

ipcopy!(dest, src) = copy!(dest, src)
ipcopy!(dest::SVector, src) = src

function qfit_min(xfm, xf0, xfp)
    xm, fm = xfm
    x0, f0 = xf0
    xp, fp = xfp
    @assert(xp > x0 && x0 > xm && isfinite(xm) && isfinite(xp))
    cm = fm/((xm-x0)*(xm-xp))  # coefficients of Lagrange polynomial
    c0 = f0/((x0-xm)*(x0-xp))
    cp = fp/((xp-xm)*(xp-x0))
    csum = cm+c0+cp
    if csum > 0
        # It's convex
        return (cm*(x0+xp) + c0*(xm+xp) + cp*(xm+x0))/(2*csum)
    end
    if fm == f0 == fp
        return x0
    end
    @assert(f0 >= fm || f0 >= fp)
    # This should only happen at a box edge. Just return something
    # to indicate the direction of the minimum
    if f0 > fm && fp > fm
        return xm - (x0 - xm)
    end
    return xp + (xp - x0)
end
