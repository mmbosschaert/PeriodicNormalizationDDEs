module DDEBifTool

import Base.vec
import Base.*
import Base.+
import Base.-

# overload the NaNMath.cos function to allow for symbolic Num
# used in build_function in  characteristic_matrix_unevaluated.jl
export characteristic_matrices, getJet, eigenfunctionsBT, BTnormalformcoefficients
export coefficient_matrices, characteristic_matrices_unevaluated
export characteristic_matrices_unevaluated_re_im
export Hopf_res!
export Hopf, stability
export newton
export stst
export vec
export locate_double_hopf
export locate_zero_hopf
export detect_genh
export normalform
export psol
export PsolLPC
export psol_res
export dL
export -, *
export defsystem_psol_jac
export defsystem_fold_jac_β
export psolLPC_res
export psolLPC_res_β
export interpolate
export d_interpolate
export fold_q1_approx
export fold_tangent
export remesh_psol
export defsystem_psol_jac_no_mod
export multipliers
export psol_pd_res
export pd_q1_approx
export defsystem_pd_jac
export psol_pd
export psol_ns
export psol_ns_res
export psol_to_ns
export differential_equation_part_complex
export ns_w_approx
export ns_res_jac
export psol_ns_res_broatcast
export pd_w_approx
export psol_to_pd
export psol_pd_res_as_ns
export pd_res_jac
export continue!
export SetupStstBranch
export LocateHopfPoints
export SetupHopfBranch
export reverse_branch!
export SetupPsolBranch
export multipliers!
export detectSpecialPoints
export SetupLPCBranch
export normal_form_coefficients!
export SetupNSBranch
export nmfm_coefficient
export SetupPDBranch
export LocateGPDPoints
export point_to_hoho
export hopf_from_hoho
export get_params
export get_nmfm_coefficients
export point_to_genhopf
export locate_genh
export +, *
export doubleHopfToPsol
export generalizedHopfToPsol
export LocateDoubleHopf

using LinearAlgebra
using Symbolics
using Polynomials
using NonlinearEigenproblems
using GaussQuadrature
using NLsolve
using Setfield
using SparseArrays
using NaNMath
# NaNMath.cos(x::Num) = cos(x)
# NaNMath.sin(x::Num) = sin(x)

include("./characteristic_matrix.jl")
include("./characteristic_matrix_unevaluated.jl")
include("./coefficient_matrices.jl")
include("./getJet.jl")
include("./eigenfunctionsBT.jl")
include("./BTnormalformcoefficients.jl")
include("./hopf.jl")
include("./newton.jl")
include("./continuation.jl")
include("./stst.jl")
include("./branch.jl")
include("./stability.jl")
include("./doubleHopf.jl")
include("./zeroHopf.jl")
include("./psol.jl")
include("./lpc.jl")
include("./ns.jl")
include("./multipliers.jl")
include("./interpolation.jl")
include("./pd.jl")
include("./borderedInverse.jl")
include("./coefficients.jl")
include("./multilinear_forms.jl")
include("./NewtonCotesWeights.jl")
include("./generalizedHopf.jl")

# Auxiliary functions
# auxiliary function to extract parameters
get_params(points) = hcat([p.parameters for p in points]...)
# auxiliary function to extract normal form coefficients
get_nmfm_coefficients(branch) = [p.nmfm for p in branch.points]

end # module
