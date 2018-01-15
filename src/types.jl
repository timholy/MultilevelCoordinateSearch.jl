"""
    root = Box{T}()
    child = Box(parent, level, edge)

Create a new box, a structure describing a rectangular domain nested
within a larger region. `parent` is the parent `Box`, `level` is an
integer assigned to indicate the likely priority that this box should be split,
and `edge` stores the value of one particular coordinate defining
the boundary between this box and its "partner" as defined by a golden
section split (or the edge of the outer bounds if this box lies between an
evaluation point and the outer bounds).

The clearest way to understand the design of `Box`
is by reference to Fig 2 in the [mcs](@ref) citation.
"""
mutable struct Box{T}
    level::Int
    edge::T  # the edge not corresponding to one of the parent's xvalues (splits that occur between dots in Fig 2)
    parent_cindex::Int # of its parent's children, which one is this?
    splitdim::Int # the dimension along which this box has been split (0 for leaf-boxes)
    parent::Box{T}
    xvalues::Vector{T} # the values of x_splitdim at which f is evaluated
    fvalues::Vector{T} # the corresponding values of f
    children::Vector{Box{T}}

    function Box{T}() where T
        # Create the root box
        box = new{T}(1, zero(T), typemax(Int)-1, 0)  # the root box will always be split along dimension 1
        box.parent = box
        return box
    end
    function Box{T}(parent::Box, level::Int, edge::T) where T
        # Create a new child and store it in the parent
        if !isdefined(parent, :children)
            parent.children = Box{T}[]
        end
        parent_cindex = length(parent.children) + 1
        box = new{T}(level, edge, parent_cindex, 0, parent)
        push!(parent.children, box)
        box
    end
end

Box(parent::Box{T}, splitdim, level, edge) where T =
    Box{T}(parent, splitdim, level, edge)

isroot(box::Box) = box.parent == box
issplit(box::Box) = isdefined(box, :children)
