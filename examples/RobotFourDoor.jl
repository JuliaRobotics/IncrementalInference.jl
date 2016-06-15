using IncrementalInference, KernelDensityEstimate
using Gadfly # for vstack
using Cairo, Fontconfig # for drawing PNG/PDF

# if you want to add multiple process here rather than calling $ julia -p4
# addprocs(3)

gt = Dict{ASCIIString, Array{Float64,2}}()
# HMM computed ground truth for first 3 poses only
gt["x1"]=[[-100.0; 1.96]';[0.0; 1.96]']'
gt["x2"]=[[-50.0; 3.1]';[50.0; 3.1]']'
gt["x3"]=[[100.0; 3.05]';[0.0; 3.05]']'

fg = emptyFactorGraph()

N=100

doors = [-100.0;0.0;100.0;300.0]'
pd = kde!(doors,[3.0])
pd = resample(pd,N);
bws = getBW(pd)[:,1]
doors2 = getPoints(pd);
v1 = addNode!(fg,"x1",doors,N=N)
f1  = addFactor!(fg,Obsv2([v1], doors2, bws', [1.0]))

tem = 2.0*randn(1,N)+getVal(v1)+50.0
v2 = addNode!(fg,"x2", tem, N=N)
addFactor!(fg,Odo([v1;v2],[50.0]',[2.0]',[1.0]))

# monocular sighting would look something like
#addFactor!(fg, Mono, ["x3","l1"], [14.0], [1.0], [1.0])
#addFactor!(fg, Mono, ["x4","l1"], [11.0], [1.0], [1.0])

v3=addNode!(fg,"x3",4.0*randn(1,N)+getVal(v2)+50.0, N=N)
addFactor!(fg,Odo([v2;v3],[50.0]',[4.0]',[1.0]))
f2 = addFactor!(fg,Obsv2([v3], doors2, bws', [1.0]))

if true

    v4=addNode!(fg,"x4",2.0*randn(1,N)+getVal(v3)+50.0, N=N)
    addFactor!(fg,Odo([v3;v4],[50.0]',[2.0]',[1.0]))

    if true
        l1=addNode!(fg, "l1", 0.5*randn(1,N)+getVal(v3)+64.0, N=N)
        addFactor!(fg, Ranged([v3,l1], [64.0],[0.5],[1.0]))
        addFactor!(fg, Ranged([v4,l1], [16.0],[0.5],[1.0]))
    end


    v5=addNode!(fg,"x5",2.0*randn(1,N)+getVal(v4)+50.0, N=N)
    addFactor!(fg,Odo([v4;v5],[50.0]',[2.0]',[1.0]))


    v6=addNode!(fg,"x6",1.25*randn(1,N)+getVal(v5)+40.0, N=N)
    addFactor!(fg,Odo([v5;v6],[40.0]',[1.25]',[1.0]))


    v7=addNode!(fg,"x7",2.0*randn(1,N)+getVal(v6) +60.0, N=N)
    addFactor!(fg,Odo([v6;v7],[60.0]',[2.0]',[1.0]))

    f3 = addFactor!(fg,Obsv2([v7], doors, bws', [1.0]))


    # HMM computed ground truth, extended for 7 poses with landmark
    gt = Dict{ASCIIString, Array{Float64,2}}()
    gt["x1"]=([0.0;1.97304 ]')' # -0.0342366
    gt["x2"]=([50.0; 2.83153 ]')' # 49.8797
    gt["x3"]=([100.0; 1.65557 ]')' # 99.8351
    gt["x4"]=([150.0; 1.64945 ]')' # 148.637
    gt["x5"]=([200.0; 1.77992 ]')' # 198.62
    gt["x6"]=([240.0; 2.20466 ]')' # 238.492
    gt["x7"]=([300.0; 2.14353 ]')' # 298.467
    gt["l1"]=([165.0; 1.17284 ]')' # 164.102
end

writeGraphPdf(fg);
tree = prepBatchTree!(fg,drawpdf=true);


# list vertices in fg
@show xx,ll = ls(fg)

# do belief propagation inference over tree once
[inferOverTree!(fg, tree) for i in 1:1];

# draw all beliefs
DOYTICKS = false
xx,ll = ls(fg)
msgPlots = drawHorBeliefsList(fg, xx, gt=gt,nhor=2);
evalstr = ""
for i in 1:length(msgPlots)
    evalstr = string(evalstr, ",msgPlots[$(i)]")
end
pl = eval(parse(string("vstack(",evalstr[2:end],")")));
Gadfly.draw(PDF("4doors.pdf",15cm,20cm),pl) # can also do PNG
# Gadfly.set_default_plot_size(15cm, 20cm)
# pl


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



using ProtoBuf

iob = PipeBuffer()
pd = convert(PackedVariableNodeData,fg.v[1].attributes["data"])
writeproto(iob, pd)
dd = readproto(iob, PackedVariableNodeData());
unpckd = convert(VariableNodeData, dd)
