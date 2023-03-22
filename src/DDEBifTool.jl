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
include("./multipliers.jl")
include("./interpolation.jl")

end # module
