
function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
	return (nlsolve(   (X, res) -> residFnc(res, measurement, parameters, X), x0 )).zero
end


function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
		# n = length(X)
    Y = deepcopy(Alt)
    for i in 1:d
			Y[p[i]] = X[i]
		end
    return Y
end

# residual function must have the form residFnc!(res, x)
function numericRootGenericRandomizedFncOld(
      residFnc!::Function,
      zDim::Int,
      xDim::Int,
      x0::Vector{Float64};
      perturb::Float64=0.00001  )

	# xDim = length(x0)
	if zDim < xDim
		p = collect(1:xDim);
		shuffle!(p);
		# p1 = p.==1; p2 = p.==2; p3 = p.==3
    # TODO -- refactor to use functors instead
		r = nlsolve(    (x, res) -> residFnc!(shuffleXAltD(x, x0, zDim, p), res ),
                    x0[p[1:zDim]] + perturb*randn(zDim)
               )

    return shuffleXAltD(r.zero, x0, zDim, p );
	else
    # return (nlsolve(  (x, res) -> residFnc!(res, x),
    #                   x0      )).zero
    # @show x0
    return ( nlsolve(  residFnc!, x0 + perturb*randn(zDim) ) ).zero
	end
end

type FastGenericRoot <: Function
  p::Vector{Int}
  perturb::Vector{Float64}
  x0::Vector{Float64}
  Y::Vector{Float64}
  xDim::Int
  zDim::Int
  usrfnc::Function
  FastGenericRoot(xDim::Int, zDim::Int, residfnc::Function) =
      new(collect(1:xDim), zeros(xDim), zeros(xDim), zeros(xDim), xDim, zDim, residfnc)
end

function shuffleXAltD!(fr::FastGenericRoot, X::Vector{Float64})
  copy!(fr.Y, fr.x0)
  for i in 1:fr.zDim
    fr.Y[fr.p[i]] = X[i]
  end
  # @show fr.Y
  nothing
end

function (fr::FastGenericRoot)( x::Vector{Float64}, res::Vector{Float64} )
  #
  shuffleXAltD!(fr, x)  #(fr.Y, x, fr.x0, fr.zDim, fr.p)
  fr.usrfnc( fr.Y, res)
end

function numericRootGenericRandomizedFnc(
      residFnc!::Function,
      zDim::Int,
      xDim::Int,
      x0::Vector{Float64};
      perturb::Float64=0.00001,
      testshuffle::Bool=false   )


	# xDim = length(x0)
	if zDim < xDim || testshuffle
    # TODO -- this only start of refactoring for inplace, more to come
    wrappedresid = FastGenericRoot(xDim, zDim, residFnc!)
    shuffle!(wrappedresid.p);
    wrappedresid.perturb = perturb*randn(zDim)
    wrappedresid.x0 = x0
    r = nlsolve(    wrappedresid,
                    x0[wrappedresid.p[1:zDim]] + wrappedresid.perturb
               )
    # r = nlsolve(    (x, res) -> residFnc!(shuffleXAltD(x, x0, zDim, p), res ),
    #                x0[wrappedresid.p[1:zDim]] + wrappedresid.perturb
    #            )

    return shuffleXAltD(r.zero, x0, zDim, wrappedresid.p );
	else
    # return (nlsolve(  (x, res) -> residFnc!(res, x),
    #                   x0      )).zero
    # @show x0
    return ( nlsolve(  residFnc!, x0 + perturb*randn(zDim) ) ).zero
	end
end






#
