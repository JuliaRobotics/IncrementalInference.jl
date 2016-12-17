
function numericRoot(residFnc::Function, measurement, parameters, x0::Vector{Float64})
	return (nlsolve(   (X, res) -> residFnc(res, measurement, parameters, X), x0 )).zero
end


function shuffleXAltD(X::Vector{Float64}, Alt::Vector{Float64}, d::Int, p::Vector{Int})
		n = length(X)
    Y = deepcopy(Alt)
		for i in 1:d
			Y[p[i]] = X[i]
		end
    return Y
end

# residual function must have the form residFnc!(res, x)
function numericRootGenericRandomizedFnc(
      residFnc!::Function,
      zDim::Int,
      xDim::Int,
      x0::Vector{Float64},
      perturb::Float64=0.001  )

	# xDim = length(x0)
	if zDim < xDim
		p = collect(1:xDim);
		shuffle!(p);
		# p1 = p.==1; p2 = p.==2; p3 = p.==3
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
