
# ==============================================================================
#  A few methods to upstream
# FIXME UPSTREAM TO DFG
# ==============================================================================

DistributedFactorGraphs.getStateKind(kind::StateType) = kind

"""
    $(SIGNATURES)

Get the state belief estimate, which is of type ::HomotopyDensity.
"""
function getBelief(state::State; newbw::Bool = true)
  state.belief
  # return manikde!(getManifold(getStateKind(vnd)), getVal(vnd); bw = getBW(vnd)[:, 1], newbw)
end
function getBelief(v::VariableCompute, solvekey::Symbol = :default; newbw::Bool = true)
  return getBelief(getState(v, solvekey); newbw)
end
function getBelief(dfg::AbstractDFG, lbl::Symbol, solvekey::Symbol = :default; newbw::Bool = true)
  return getBelief(getVariable(dfg, lbl), solvekey; newbw)
end









## DEPRECATE BELOW

# ==============================================================================
#  Topology types that specialize AbstractHomotopyTopology (defined in DFG)
# ==============================================================================
# L1 structural nodes only. No L2 samples. (Schema: `means`, `weights`, `forms` populated. `points` empty.)
struct RootsOnlyTopology <: DFG.AbstractHomotopyTopology end
# L2 raw samples only. No L1 structure. (Schema: `points`, `bandwidths` populated. `means` empty.)
struct LeavesOnlyTopology <: DFG.AbstractHomotopyTopology end

# Convenience constructors for HomotopyDensityDFG with topology dispatch
function DFG.HomotopyDensityDFG(::LeavesOnlyTopology, T::DFG.AbstractStateType; kwargs...)
    dim = DFG.getDimension(T)
    return DFG.HomotopyDensityDFG{typeof(T), DFG.getPointType(T)}(;
        reprkind = DFG.HomotopyReprDFG(LeavesOnlyTopology(), DFG.DefaultFormKind(), T, nothing),
        trailing_forms = sparsevec(Dict(1 => zeros(dim, dim))),
        kwargs...,
    )
end

function DFG.HomotopyDensityDFG(::RootsOnlyTopology, T::DFG.AbstractStateType; kwargs...)
    return DFG.HomotopyDensityDFG{typeof(T), DFG.getPointType(T)}(;
        reprkind = DFG.HomotopyReprDFG(RootsOnlyTopology(), DFG.DefaultFormKind(), T, nothing),
        kwargs...,
    )
end