function continuation!(branch, v₀, n, f, df; number_of_points = 300, δ = 1e-03, parameterbounds = nothing)

    xcorrected = vec(branch[end])

    for _ = 2:number_of_points
        xnew =  xcorrected + δ*v₀
        xcorrected, _, v₀ = newton(f, df, xnew, xcorrected, v₀, n, tol = 1e-14)

        # calulate stability
        push!(branch,vec_to_point(typeof(branch[1]),xcorrected,n))

        # check if piont is within bounds
        if parameterbounds !== nothing 
            if !all(>(0),(parameterbounds.min .< branch[end].parameters .< parameterbounds.max))
                break
            end
        end

    end

end

function continue_psol(jet,f,df,branch;δ=1,max_number_of_points=100,δmin=0.001,δmax=2)
    i = 1
    while true
        psol_guess = branch[end] + δ*(branch[end] - branch[end-1])
        x₀ = [vcat(psol_guess.profile...); psol_guess.parameters[2]; psol_guess.period ...]
        # sol = nlsolve(x->f(x,psol_guess), x₀, iterations=20)
        sol = nlsolve(x->f(x,psol_guess), x-> df(x,psol_guess), x₀, iterations=20)
        if sol.f_converged
            psol_corrected = psol([reshape(sol.zero[1:end-2],2,:)[:,j] for j ∈ 1:length(psol_guess.mesh)],[psol_guess.parameters[1]; sol.zero[end-1]], psol_guess.mesh, sol.zero[end], psol_guess.ncol, nothing, nothing)
            push!(branch, psol_corrected)
            i += 1
            δ = sol.iterations > 2 ? δ/1.1 : δ*1.1 
        else
            δ /= 1.1
        end
        if i > max_number_of_points
            break
        end
        if δ > δmax
            δ = δmax
        end
        if δ < δmin
            println("Minimum stepsize")
            break
        end
    end
    branch
end


function continue_psol_fold(jet,f,df,branch;δ=1,max_number_of_points=100,δmin=0.001,δmax=2)
    i = 1

    ncol = branch[1].ncol
    ntst = convert(Int,(length(branch[1].profile)-1)/ncol)
    dims = length(branch[1].profile[1])

    while true
        @show i
        psol_guess = branch[end] + δ*(branch[end] - branch[end-1])
        x₀ = [vcat(psol_guess.eigenvector...); vcat(psol_guess.profile...); psol_guess.parameters; psol_guess.period ...]
        sol = nlsolve(x->f(x,psol_guess), x₀, iterations=20)
        # sol = nlsolve(x->f(x,psol_guess), x-> df(x,psol_guess), x₀, iterations=20)
        if sol.f_converged
            psol_corrected = psol_fold([reshape(sol.zero[dims*(ntst*ncol+1)+1:2*dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1],
                  [reshape(sol.zero[1:dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1],
                  [sol.zero[end-2]; sol.zero[end-1]], psol_guess.mesh, sol.zero[end], 
                  psol_guess.ncol, nothing, nothing,0.0)
            push!(branch, psol_corrected)
            i += 1
            δ = sol.iterations > 2 ? δ/1.1 : δ*1.1 
        else
            δ /= 1.1
        end
        if i > max_number_of_points
            break
        end
        if δ > δmax
            δ = δmax
        end
        if δ < δmin
            println("Minimum stepsize")
            break
        end
    end
    branch
end

function continue_psol_fold_β(f,df,branch;δ=1,max_number_of_points=100,δmin=0.001,δmax=2)
    i = 1

    ncol = branch[1].ncol
    ntst = convert(Int,(length(branch[1].profile)-1)/ncol)
    dims = length(branch[1].profile[1])

    while true
        fold_guess = branch[end] + δ*(branch[end] - branch[end-1])

        x₀ = [vcat(fold_guess.profile...); vcat(fold_guess.eigenvector...); fold_guess.parameters; fold_guess.period; fold_guess.beta ...]
        sol = nlsolve(x->f(x,fold_guess),x->df(x,fold_guess), x₀, iterations=20)
        if sol.f_converged
            fold_corrected = psol_fold([reshape(sol.zero[1:dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1],
                              [reshape(sol.zero[dims*(ntst*ncol+1)+1:2dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1],
                              sol.zero[end-3:end-2], fold_guess.mesh, sol.zero[end-1], 
                              fold_guess.ncol, nothing, nothing, sol.zero[end])
            push!(branch, fold_corrected)
            i += 1
            δ = sol.iterations > 2 ? δ/1.1 : δ*1.1 
        else
            δ /= 1.1
        end
        if i > max_number_of_points
            break
        end
        if δ > δmax
            δ = δmax
        end
        if δ < δmin
            println("Minimum stepsize")
            break
        end
    end
    branch
end

function continue_psol!(branch,f,df,V;δ=1,MaxNumberofSteps=100,δmin=0.001,δmax=1.2)

    psol_guess = branch.points[end]
    ncol = psol_guess.ncol
    ntst = convert(Int,(length(psol_guess.profile)-1)/ncol)
    dims = length(psol_guess.profile[1])

    iterations=0
    for step = 1:MaxNumberofSteps

        # step != 1 || length(branch.points) != 0 ?  x₀ += δ₀*V : x₀
        psol_guess = branch.points[end]
        x₀ = [vcat(psol_guess.profile...); psol_guess.parameters[2]; psol_guess.period ...]
        x₀ += δ*V # new guess
        x₀, iterations, V, converged  = newton(f, df, x₀, V, branch.points[end]; tol=1e-8)

        if ~converged
            println("No convergence")
            break
        end


        # condition on convergence here
        iterations < 5 ? δ *= 1.2 : δ /= 1.2

        psol_corrected = psol([reshape(x₀[1:end-2],2,:)[:,j] for j in 1:length(psol_guess.mesh)], 
                              [psol_guess.parameters[1]; x₀[end-1]], 
                              psol_guess.mesh, x₀[end], ncol, nothing, nothing)

        push!(branch.points, psol_corrected)
        push!(branch.tangents, V)
        push!(branch.stepsizes, δ)

        # if step % 25 == 0
        #     println("Generation new mesh")
        #     τ, hom, vpsnew = adaptmeshHom(reshape(x₀[1:dims*(N*m+1)], dims, N*m+1),
        #                                   reshape(V[1:dims*(N*m+1)], dims, N*m+1), τ, N, m)

        #     V[1:dims*(N*m+1)] = vpsnew

        #     # update solution
        #     x₀[1:dims*(N*m+1)] = hom
        #     branch.points[end] = x₀
        #     branch.meshes[end] = τ
        # end

        println("Step $step: iterations: $iterations stepsize: $δ")

        if δ < 1e-15
            println("Stepsize to small: δ = $δ")
            break
        end
    end
end

function continue_psol_fold_β!(branch,f,df,V;δ=1,MaxNumberofSteps=100,δmin=0.001,δmax=1.2)

    psol_guess = branch.points[end]
    ncol = psol_guess.ncol
    ntst = convert(Int,(length(psol_guess.profile)-1)/ncol)
    dims = length(psol_guess.profile[1])

    iterations=0
    for step = 1:MaxNumberofSteps

        # step != 1 || length(branch.points) != 0 ?  x₀ += δ₀*V : x₀
        psol_guess = branch.points[end]
        x₀ = [vcat(psol_guess.profile...); vcat(psol_guess.eigenvector...); psol_guess.parameters; psol_guess.period; psol_guess.beta ...]
        x₀ += δ*V # new guess
        x₀, iterations, V, converged  = newton(f, df, x₀, V, branch.points[end]; tol=1e-8)

        if ~converged
            println("No convergence")
            break
        end


        # condition on convergence here
        iterations < 5 ? δ *= 1.2 : δ /= 1.2

        fold_corrected = psol_fold([reshape(x₀[1:dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1],
                          [reshape(x₀[dims*(ntst*ncol+1)+1:2dims*(ntst*ncol+1)],2,:)[:,j] for j in 1:ntst*ncol+1],
                          x₀[end-3:end-2], branch.points[end].mesh, x₀[end-1], 
                          branch.points[end].ncol, nothing, nothing, x₀[end])

        push!(branch.points, fold_corrected)
        push!(branch.tangents, V)
        push!(branch.stepsizes, δ)

        # if step % 25 == 0
        #     println("Generation new mesh")
        #     τ, hom, vpsnew = adaptmeshHom(reshape(x₀[1:dims*(N*m+1)], dims, N*m+1),
        #                                   reshape(V[1:dims*(N*m+1)], dims, N*m+1), τ, N, m)

        #     V[1:dims*(N*m+1)] = vpsnew

        #     # update solution
        #     x₀[1:dims*(N*m+1)] = hom
        #     branch.points[end] = x₀
        #     branch.meshes[end] = τ
        # end

        println("Step $step: iterations: $iterations stepsize: $δ")

        if δ < 1e-15
            println("Stepsize to small: δ = $δ")
            break
        end
    end
end


