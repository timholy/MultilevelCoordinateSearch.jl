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

@testset "Init: box edges extend to region edges" begin
    xstar, fstar, box = MultilevelCoordinateSearch.init_boxes(camel, [0.0,0.0], ([-3.0,0.0,3.0], [-2.0,0.0,2.0]), u, v)
    @test xstar == [0.0, 0.0]
    @test fstar == camel(xstar)
    @test length(box.parent.children) == 4 && length(box.parent.parent.children) == 4
    r = MultilevelCoordinateSearch.get_root(box)
    @test r.fvalues[2] < min(r.fvalues[1], r.fvalues[3]) && [r.children[i].level for i = 1:4] == [3, 2, 2, 3]
    b = box.parent
    @test b.fvalues[2] < min(b.fvalues[1], b.fvalues[3]) && [b.children[i].level for i = 1:4] == [4, 3, 3, 4]
    @test MultilevelCoordinateSearch.splitindexes(b.children[1]) == (1, 2, 1)
    @test MultilevelCoordinateSearch.splitindexes(b.children[2]) == (1, 2, 2)
    @test MultilevelCoordinateSearch.splitindexes(b.children[3]) == (2, 3, 2)
    @test MultilevelCoordinateSearch.splitindexes(b.children[4]) == (2, 3, 3)
    pqs = MultilevelCoordinateSearch.init_sweep_queues(box, 4)
    @test isempty(pqs[1])
    @test length(pqs[2]) == 1
    @test length(pqs[3]) == 4
    @test length(pqs[4]) == 2
end

@testset "Init: box edges interior to region edges" begin
    xstar, fstar, box = MultilevelCoordinateSearch.init_boxes(camel, [0.0,0.0], ([-2.0,0.0,2.0], [-1.0,0.0,1.0]), u, v)
    @test xstar == [0.0, 0.0]
    @test fstar == camel(xstar)
    @test length(box.parent.children) == 6 && length(box.parent.parent.children) == 6
    r = MultilevelCoordinateSearch.get_root(box)
    @test r.fvalues[2] < min(r.fvalues[1], r.fvalues[3]) && [r.children[i].level for i = 1:6] == [2, 3, 2, 2, 3, 2]
    b = box.parent
    @test b.fvalues[2] == min(b.fvalues[1], b.fvalues[3]) && [b.children[i].level for i = 1:6] == [3, 4, 3, 4, 3, 3]
    @test MultilevelCoordinateSearch.splitindexes(b.children[1]) == (0, 1, 1)
    @test MultilevelCoordinateSearch.splitindexes(b.children[2]) == (1, 2, 1)
    @test MultilevelCoordinateSearch.splitindexes(b.children[3]) == (1, 2, 2)
    @test MultilevelCoordinateSearch.splitindexes(b.children[4]) == (2, 3, 2)
    @test MultilevelCoordinateSearch.splitindexes(b.children[5]) == (2, 3, 3)
    @test MultilevelCoordinateSearch.splitindexes(b.children[6]) == (3, 4, 3)
    pqs = MultilevelCoordinateSearch.init_sweep_queues(box, 4)
    @test isempty(pqs[1])
    @test length(pqs[2]) == 3
    @test length(pqs[3]) == 6
    @test length(pqs[4]) == 2
end

@testset "Init: box edges half to region edges" begin
    xstar, fstar, box = MultilevelCoordinateSearch.init_boxes(camel, [0.0,0.0], ([-3.0,0.0,2.0], [-1.0,0.0,2.0]), u, v)
    @test xstar == [0.0, 0.0]
    @test fstar == camel(xstar)
    @test length(box.parent.children) == 5 && length(box.parent.parent.children) == 5
    r = MultilevelCoordinateSearch.get_root(box)
    @test r.fvalues[2] < min(r.fvalues[1], r.fvalues[3]) && [r.children[i].level for i = 1:5] == [3, 2, 2, 3, 2]
    b = box.parent
    @test (b.fvalues[2] == b.fvalues[1] && b.fvalues[2] < b.fvalues[3]) && [b.children[i].level for i = 1:5] == [3, 4, 3, 3, 4]
    @test MultilevelCoordinateSearch.splitindexes(r.children[1]) == (1, 2, 1)
    @test MultilevelCoordinateSearch.splitindexes(r.children[2]) == (1, 2, 2)
    @test MultilevelCoordinateSearch.splitindexes(r.children[3]) == (2, 3, 2)
    @test MultilevelCoordinateSearch.splitindexes(r.children[4]) == (2, 3, 3)
    @test MultilevelCoordinateSearch.splitindexes(r.children[5]) == (3, 4, 3)
    @test MultilevelCoordinateSearch.splitindexes(b.children[1]) == (0, 1, 1)
    @test MultilevelCoordinateSearch.splitindexes(b.children[2]) == (1, 2, 1)
    @test MultilevelCoordinateSearch.splitindexes(b.children[3]) == (1, 2, 2)
    @test MultilevelCoordinateSearch.splitindexes(b.children[4]) == (2, 3, 2)
    @test MultilevelCoordinateSearch.splitindexes(b.children[5]) == (2, 3, 3)
    pqs = MultilevelCoordinateSearch.init_sweep_queues(box, 4)
    @test isempty(pqs[1])
    @test length(pqs[2]) == 2
    @test length(pqs[3]) == 5
    @test length(pqs[4]) == 2
end
