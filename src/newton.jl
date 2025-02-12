# function newton(jet, f, df, x₀, V, τ, dims, x₀prev, Qprev, complexIndx; tol=1e-8)
# 
#     maxIter = 100
# 
#     iterations = 0
#     for i=1:maxIter
# 
#         hompred0, α0, eq, Q, λ = rearrangesolutionFullSystem(x₀, N, m, dims)
#     
#             
#         QQ = [BVPHom(jet, hompred0, eq, α0, λ, Q, homPrevious, Qprev, N, m, T, τ, complexIndx); 0]
# 
#         # @show norm(QQ)
#         if norm(QQ) < tol
#             iterations = i
#             break
#         end
#         F = zeros(dims*(N*m+1)+1+dims + dims*(dims+1) + 1, dims*(N*m+1)+2+dims + dims*(dims+1))
#         BVPHomJac!(F, jet, hompred0, eq, α0, λ, Q, homPrevious, Qprev, N, m, T, τ, complexIndx)
#         F[end,:] = V'
# 
#         D = sparse(F)\[QQ [zeros(length(V)-1); 1]]
#         dx = D[:,1]
#         V = D[:,2]/norm(D[:,2])
#         x₀ = x₀ - dx
# 
#     end
# 
#     x₀, iterations, V
# end

function newton(f, df, x₀, x₀prev, v₀, n; tol=1e-10, maxIter=100)
  iterations = 0
  for i = 1:maxIter
    fx = f(x₀, x₀prev)
    if norm(fx) < tol
      iterations = i
      break
    end

    jac = df(x₀, x₀prev)
    D = [jac; v₀'] \ [[fx; 0.0] [zeros(3n + 2); 1.0]]
    dx = D[:, 1]
    v₀ = D[:, 2] / norm(D[:, 2])
    x₀ = x₀ - dx
  end
  x₀, iterations, v₀
end

function newton(f, df, x₀, x₀prev; tol=1e-10, maxIter=100)
  iterations = 0
  converged = false
  for i = 1:maxIter
    fx = f(x₀, x₀prev)
    if norm(fx) < tol
      iterations = i
      converged = true
      break
    end

    jac = df(x₀, x₀prev)
    dx = jac \ fx
    x₀ = x₀ - dx
  end
  x₀, iterations, converged
end

function newton(f, df, x₀, V, x_ref; tol=1e-8, maxiterations=20)
  iterations = 0
  converged = false
  for i = 1:maxiterations
    fx = f(x₀, x_ref)
    # @show norm(fx)
    if norm(fx) < tol
      iterations = i
      converged = true
      break
    end

    jac = df(x₀, x_ref)
    #jac_df = Main.FiniteDiff.finite_difference_jacobian(x -> f(x, x_ref), x₀)

    # @show findall(abs.(jac - jac_df) .>= 1e-04)

    jac[end, :] = V'

    if !all(isfinite, jac) || !all(isfinite, fx) || !all(isfinite, V)
        println("jac, fx, or V contains NaN/Inf; cannot proceed with solve.")
        return x₀, iterations, V, converged
    end

    # check if jac is non-singular
    if det(jac) ≈ 0.0
        println("Jacobian is singular; cannot proceed with solve.")
        return x₀, iterations, V, converged
    end

    D = jac \ [fx [zeros(length(V) - 1); 1.0]]
    dx = D[:, 1]
    V = D[:, 2] / norm(D[:, 2])
    x₀ = x₀ - dx
  end

  x₀, iterations, V, converged
end
