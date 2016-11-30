# test Pose2DPoint2D constraint evaluation function

using IncrementalInference
using Base.Test

begin
N = 100
fg = emptyFactorGraph()
fg.sessionname="NA"


initCov = diagm([0.03;0.03;0.001])
odoCov = diagm([3.0;3.0;0.01])

# Some starting position
v1 = addNode!(fg, "x1", zeros(3,1), diagm([1.0;1.0;0.1]), N=N)

# Two landmarks
l1 = addNode!(fg, "l1", ([10.0;0.0]')', diagm([1.0;1.0]), N=N)

# and constraints to pose x1
rhoZ1 = norm([10.0;0.0])
ppr = Pose2DPoint2DRange([rhoZ1], 2.0, [1.0])

f1 = addFactor!(fg, [v1;l1], ppr)


pts = evalFactor2(fg, f1, l1.index)
# @show sum(sqrt(sum(pts.^2, 1 )) .< 10.0)
# @test sum(sqrt(sum(pts.^2, 1 )) .< 10.0) == 0

v2 = addNode!(fg, "x2", ([50.0;0.0;0.0]')', diagm([1.0;1.0;0.05]), N=N)
ppc = Pose2Pose2(([50.0;0.0;0.0]')', 0.1*eye(3), [1.0])
f2 = addFactor!(fg, [v1;v2], ppr)
pts = evalFactor2(fg, f2, v1.index)

@test norm(Base.mean(pts,2)[1:2]-[50.0;0.0]) < 10.0
@test abs(Base.mean(pts,2)[3]) < 0.5


# add a prior somewhere
pp2 = PriorPoint2D([10.0;0.0], diagm([1.0;1.0]), [1.0])

f3 = addFactor!(fg,[l1], pp2)
pts = evalFactor2(fg, f3, l1.index)

@test norm(Base.mean(pts,2)[:]-[10.0;0.0]) < 3.0
end
