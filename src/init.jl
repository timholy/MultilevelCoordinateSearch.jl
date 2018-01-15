function init_boxes(f, x0::AbstractVector{T}, splits, u::AbstractVector{T}, v::AbstractVector{T}) where T
    xstar = copy(x0)
    fstar = f(xstar)
    root = box = Box{T}()
    n = length(x0)
    length(splits) == n || throw(DimensionMismatch("need one vector of splits per dimension, got $(length(splits)) for $n dimensions"))
    xtmp = copy(xstar)
    for i = 1:n
        xtmp = ipcopy!(xtmp, xstar)
        xsplit = splits[i]
        L = length(xsplit)
        L >= 3 || error("there have to be at least 3 evaluations per coordinate")
        issorted(xsplit) || error("splits must be in increasing order")
        # Evaluate f along the splits, keeping track of the best
        fsplit = Vector{typeof(fstar)}(uninitialized, L)
        fmin, idxmin = oftype(fstar, Inf), 0
        for l = 1:L
            xtmp = replacecoordinate!(xtmp, i, xsplit[l])
            ftmp = f(xtmp)
            if ftmp < fmin
                fmin = ftmp
                idxmin = l
            end
            fsplit[l] = ftmp
        end
        box.splitdim = i
        box.xvalues = copy(xsplit)
        box.fvalues = fsplit
        idxmin == 0 && error("function was not finite at any evaluation point")
        # Create the child boxes
        overhang_left = xsplit[1] > u[i]
        idx = Int(!overhang_left)
        level = box.level + 1
        while idx <= L
            if idx == 0
                # beyond the edge of xsplit
                Box{T}(box, level, u[i])  # we can drop this because it gets stored in the parent
            elseif idx == L
                if xsplit[end] < v[i]
                    Box{T}(box, level, v[i]) # also beyond edge
                end
            else
                if fsplit[idx] < fsplit[idx+1]
                    # Bigger box (with lower level) is closer to the smaller f value
                    edge = xsplit[idx] + (xsplit[idx+1]-xsplit[idx]) * ((sqrt(5)-1)/2)
                    Box{T}(box, level, edge)
                    Box{T}(box, level+1, edge)
                else
                    edge = xsplit[idx+1] - (xsplit[idx+1]-xsplit[idx]) * ((sqrt(5)-1)/2)
                    Box{T}(box, level+1, edge)
                    Box{T}(box, level, edge)
                end
            end
            idx += 1
        end
        # Update xstar
        if fmin < fstar
            fstar = fmin
            xstar = replacecoordinate!(xstar, i, xsplit[idxmin])
        end
        # Pick the best box for splitting for the next coordinate
        im, i0, ip = idxmin == 1 ? (1, 2, 3) :
                     idxmin == L ? (L-2, L-1, L) :
                     (idxmin-1, idxmin, idxmin+1)
        xqmin, _ = qfit_min((xsplit[im], fsplit[im]), (xsplit[i0], fsplit[i0]), (xsplit[ip], fsplit[ip]))
        cindex = xqmin < xsplit[im] ? (assert(idxmin == 1); 1) : # at left edge
                 xqmin > xsplit[ip] ? (assert(idxmin == L); length(box.children)) :  # at right edge
                 xqmin < xsplit[i0] ? 2*idxmin-2+overhang_left : 2*idxmin-1+overhang_left
        box = box.children[cindex]
    end
    xstar, fstar, box
end

function init_sweep_queues(lastbox::Box{T}, smax) where T
    root = get_root(lastbox)
    pqs = [PriorityQueue{Box{T},T}() for s = 1:smax]
    for box in visit_leaves(root)
        pqs[box.level][box] = minimum(box)
    end
    pqs
end

function run_sweep(pqs)
    for pq in pqs
        isempty(pq) && continue
        box = dequeue!(pq)
        if !split!(box)
            if box.level < length(pqs)
                box.level += 1
                pqs[box.level][box] = minimum(box)
            end
        end
    end
end
