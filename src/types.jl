mutable struct Box{T}
    splitdim::Int
    level::Int
    parent_cindex::Int # of its parent's children, which one is this?
    parent::Box{T}
    xvalues::Vector{T} # the values of x_splitdim at which f is evaluated
    fvalues::Vector{T} # the corresponding values of f
    children::Vector{Box{T}}

    function Box{T}() where T
        # Create the root box
        box = new{T}(1, 1, 1)  # the root box will always be split along dimension 1
        box.parent = box
        return box
    end
    function Box{T}(parent::Box, splitdim::Int, level::Int) where T
        # Create a new child and store it in the parent
        if !isdefined(parent, :children)
            parent.children = Box{T}[]
        end
        parent_cindex = length(parent.children) + 1
        box = new{T}(splitdim, level, parent_cindex, parent)
        push!(parent.children, box)
        box
    end
    # function Box{T}(level::Integer, splitdim::Integer, parent::Box, parent_cindex::Integer, xvalues, fvalues) where T
    #     return new{T}(level, Box{T}[], splitdim, xvalues, fvalues, parent, parent_cindex)
    # end
    # function Box{T}(level::Integer, splitdim::Integer, xvalues, fvalues)
    #     level == 1 || error("root box must be level 1")
    #     box = new{T}(level, Box{T}[], splitdim, xvalues, fvalues)
    #     box.parent = box
    #     box.parent_cindex = 1
    #     return box
    # end
end

# Box(level::Integer, splitdim::Integer, xvalues::Vector{T}, fvalues::Vector{T}) where T =
#     Box{T}(level, splitdim, xvalues, fvalues)
# Box(level::Integer, splitdim::Integer, xvalues::Vector{T}, fvalues::Vector{T}, parent::Box, parent_cindex::Integer) where T =
#     Box{T}(level, splitdim, xvalues, fvalues, parent, parent_cindex)

isroot(box::Box) = box.parent == box
issplit(box::Box) = isdefined(box, :children)
