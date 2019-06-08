
addprocs(3) # Remove addprocs if you want to use a single process only
using IncrementalInference
using KernelDensityEstimate, Gadfly # for vstack
using Cairo, Fontconfig # for drawing PNG/PDF


gt = Dict{Symbol, Array{Float64,2}}()
# HMM computed ground truth for first 3 poses only
gt[:x1]=[[-100.0; 1.96]';[0.0; 1.96]']'
gt[:x2]=[[-50.0; 3.1]';[50.0; 3.1]']'
gt[:x3]=[[100.0; 3.05]';[0.0; 3.05]']'

fg = initfg()

N=100

doors = reshape([-100.0;0.0;100.0;300.0],1,4)
cov = 3.0*ones(1,1)
# pd = kde!(doors,cov)
# pd = resample(pd,N);
# bws = getBW(pd)[:,1]
# doors2 = getPoints(pd);


v1 = addVariable!(fg,:x1,ContinuousScalar,N=N)
f1  = addFactor!(fg,[v1], Obsv2(doors, cov, [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addVariable!(fg,:x2, ContinuousScalar, N=N)
addFactor!(fg,[v1;v2],Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0]))

# monocular sighting would look something like
#addFactor!(fg, Mono, [:x3,:l1], [14.0], [1.0], [1.0])
#addFactor!(fg, Mono, [:x4,:l1], [11.0], [1.0], [1.0])

v3=addVariable!(fg,:x3,ContinuousScalar, N=N)
addFactor!(fg,[v2;v3],Odo(50.0*ones(1,1),4.0*ones(1,1),[1.0]))
f2 = addFactor!(fg,[v3], Obsv2(doors, cov', [1.0]))

if true

    v4=addVariable!(fg,:x4,ContinuousScalar, N=N)
    addFactor!(fg,[v3;v4],Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0]))

    if true
        l1=addVariable!(fg, :l1, ContinuousScalar, N=N)
        addFactor!(fg, [v3,l1], Ranged([64.0],[0.5],[1.0]))
        addFactor!(fg, [v4,l1], Ranged([16.0],[0.5],[1.0]))
    end


    v5=addVariable!(fg,:x5,ContinuousScalar, N=N)
    addFactor!(fg,[v4;v5],Odo(50.0*ones(1,1),2.0*ones(1,1),[1.0]))


    v6=addVariable!(fg,:x6,ContinuousScalar, N=N)
    addFactor!(fg,[v5;v6],Odo(40.0*ones(1,1),1.25*ones(1,1),[1.0]))


    v7=addVariable!(fg,:x7,ContinuousScalar, N=N)
    addFactor!(fg,[v6;v7],Odo(60.0*ones(1,1),2.0*ones(1,1),[1.0]))

    f3 = addFactor!(fg,[v7], Obsv2(doors, cov, [1.0]))


    # HMM computed ground truth, extended for 7 poses with landmark
    gt = Dict{String, Array{Float64,2}}()
    gt[:x1]=reshape(Float64[0.0;1.97304 ],2,1) # -0.0342366
    gt[:x2]=reshape(Float64[50.0; 2.83153 ],2,1) # 49.8797
    gt[:x3]=reshape(Float64[100.0; 1.65557 ],2,1) # 99.8351
    gt[:x4]=reshape(Float64[150.0; 1.64945 ],2,1) # 148.637
    gt[:x5]=reshape(Float64[200.0; 1.77992 ],2,1) # 198.62
    gt[:x6]=reshape(Float64[240.0; 2.20466 ],2,1) # 238.492
    gt[:x7]=reshape(Float64[300.0; 2.14353 ],2,1) # 298.467
    gt[:l1]=reshape(Float64[165.0; 1.17284 ],2,1) # 164.102
end

# writeGraphPdf(fg);
tree = prepBatchTree!(fg,drawpdf=true);

# spyCliqMat(tree.cliques[1])

# list vertices in fg
@show xx,ll = ls(fg)

# from = fg.v[fg.IDs[:l1]]
# to = fg.v[fg.IDs[:x7]]
# fgs = subgraphShortestPath(fg, from=from, to=to)

# using Graphs
# el = shortest_path(fg.g, ones(19),from, to)

# do belief propagation inference over tree once
[inferOverTree!(fg, tree, N=100) for i in 1:1];

# draw all beliefs
# TODO -- Update this part of the example
using RoMEPlotting

DOYTICKS = false
xx,ll = ls(fg)
msgPlots = drawHorBeliefsList(fg, xx, gt=gt,nhor=2);
evalstr = ""
for i in 1:length(msgPlots)
    evalstr = string(evalstr, ",msgPlots[$(i)]")
end
pl = eval(parse(string("vstack(",evalstr[2:end],")")));
Gadfly.draw(PNG("4doors.png",15cm,20cm),pl) # can also do PNG



# if false
#
#   # draw upward messages
#   msgPlots=drawTreeUpwardMsgs(fg, tree, N=500); # to init memory for eval(parse(string))
#   println("Upward messages for all cliques except root")
#   # vvMsgs = vstackedDensities(fg, tree, msgPlots)
#   Gadfly.set_default_plot_size(17cm, 15cm)
#   Gadfly.draw(PNG("results/testMsgs.png",17cm,15cm),vvMsgs)
#   # vvMsgs
# end
