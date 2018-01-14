using MultilevelCoordinateSearch
using Base.Test

function camel(x)
    # 6-hump camel function. Typically evaluated over [-3,3] × [-2,2].
    x1, x2 = x[1], x[2]
    x1s = x1*x1
    x2s = x2*x2
    return (4 - 2.1*x1s + x1s*x1s/3)*x1s + x1*x2 + (-4 + 4*x2s)*x2s
end

u = [-3.0,-2.0]
v = -u
splits0 = ([-2.0,0.0,2.0],
           [-1.0,0.0,1.0])
splits1 = ([-3.0,0.0,3.0],
           [-2.0,0.0,2.0])

xstar, fstar, box = MultilevelCoordinateSearch.init(camel, [0.0,0.0], splits1, u, v)
xstar, fstar, box = MultilevelCoordinateSearch.init(camel, [0.0,0.0], splits0, u, v)
