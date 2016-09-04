using Gadfly, Colors
using KernelDensityEstimate, TransformUtils, IncrementalInference

# linear product is going to use Euler angles, until improve to a manifold sampler



STD1 = Float64[]
STD2 = Float64[]
STD3 = Float64[]


SPC = [pi/4.0]
SPC = (0:62)./128.0*pi
for t in SPC
  so3sigma = 0.001
  N = 1000
  OUT = Array{Float64,2}(N,4)
  for i in 1:N
    Rp = convert(SO3, so3(so3sigma*randn(3)))
    Rx = convert(SO3, AngleAxis(t,[0;1;0]))
    aa = convert(Euler, Rx*Rp)
    # OUT[i,1] = aa.theta
    # OUT[i,2:4] = aa.ax[:]
    OUT[i,1:3] = [aa.R;aa.P;aa.Y]
  end
  if length(SPC) == 1
    plot(x=OUT[:,1], Geom.histogram)
    plot(x=OUT[:,2], Geom.histogram)
    plot(x=OUT[:,3], Geom.histogram)
  end

  s = std(OUT,1);
  # s[1] > 2*so3sigma || s[3] > 2*so3sigma ? warn("Gimbal lock becoming significant") : nothing
  push!(STD1,s[1])
  push!(STD2,s[2])
  push!(STD3,s[3])

end


plot(
layer(x=SPC*180/pi,y=STD1,Geom.line,Theme(default_color=colorant"red")),
layer(x=SPC*180/pi,y=STD2,Geom.line,Theme(default_color=colorant"black")),
layer(x=SPC*180/pi,y=STD3,Geom.line,Theme(default_color=colorant"deepskyblue")),
Coord.Cartesian(xmax=90,ymax=5*so3sigma)
)


# eval product


covOdo = 0.01*ones(6,1)
covOdo[4:6] = 0.001

mu1 = zeros(6,1)
mu1[:] = veeEuler(SE3(0))

p1 = resample(kde!(mu1, vec(covOdo)),100)

p2 = resample(kde!(mu1, vec(covOdo)),100)

p12 = p1*p2
plotKDE(marginal(p12,[6]))



#
