module DDEBifTool

import Base.vec
import Base.*
import Base.+
import Base.-

export characteristic_matrices, getJet, eigenfunctionsBT, BTnormalformcoefficients
export coefficient_matrices, characteristic_matrices_unevaluated
export characteristic_matrices_unevaluated_re_im
export Hopf_res!
export Hopf, stability
export newton
export continuation! 
export continue_psol! 
export stst
export vec
export locate_double_hopf
export locate_zero_hopf
export locate_genh
export normalform
export psol
export psol_fold
export psol_res
export dL
export -, *
export continue_psol
export defsystem_psol_jac
export defsystem_fold_jac_β
export continue_psol_fold_β!
export psol_fold_res
export psol_fold_res_β
export interpolate
export d_interpolate
export continue_psol_fold
export continue_psol_fold_β
export fold_q1_approx
export fold_tangent
export remesh_psol
export defsystem_psol_jac_no_mod
export multipliers
export psol_pd_res
export pd_q1_approx
export defsystem_pd_jac
export psol_pd
export continue_psol_pd!
export psol_ns
export psol_ns_res
export psol_to_ns
export differential_equation_part_complex
export ns_w_approx
export ns_res_jac
export psol_ns_res_broatcast
export continue_psol_ns!
export pd_w_approx
export psol_to_pd
export psol_pd_res_as_ns
export pd_res_jac

using LinearAlgebra
using Symbolics
using Polynomials
using NonlinearEigenproblems
using GaussQuadrature
using NLsolve
using Setfield
using SparseArrays

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
include("./ns.jl")
include("./multipliers.jl")
include("./interpolation.jl")
include("./pd.jl")

end # module
