function continue!(branch,f,df,V;δ=1,MaxNumberofSteps=100,δmin=0.001,δmax=1.2, parameterbounds = nothing)
    continue!(branch,f,df,vec,vec_to_point,V;δ=δ,MaxNumberofSteps=MaxNumberofSteps,δmin=δmin,δmax=δmax, parameterbounds = parameterbounds)
end

# function continue!(branch)
#     continue!(branch,branch.f,branch.df,
#         branch.point_to_vec,
#         branch.vec_to_point,
#         branch.tangents[1];
#         δ=branch.δ,
#         MaxNumberofSteps=branch.MaxNumberofSteps,
#         δmin=branch.δmin,
#         δmax=branch.δmax,
#         parameterbounds = branch.parameterbounds)
# end

function continue!(branch,f,df,point_to_vec,vec_to_point,V;δ=1,MaxNumberofSteps=250,δmin=0.001,δmax=1.2, parameterbounds = nothing, NumberOfFails = 4)

    iterations=0
    fails=0

    for step = 1:MaxNumberofSteps

        # step != 1 || length(branch.points) != 0 ?  x₀ += δ₀*V : x₀
        point_guess = branch.points[end]
        x₀ = point_to_vec(point_guess)
        x₀ += δ*V # new guess
        x₀, iterations, V, converged  = newton(f, df, x₀, V, branch.points[end]; tol=1e-8)

        if ~converged
            println("No convergence")
            fails += 1
            if fails < NumberOfFails
                δ /= 0.5
                continue
            else
                break
            end
        end
        fails = 0

        iterations < 5 ? δ *= 1.2 : δ /= 1.2
        δ = min(δ,δmax)

        pd_corrected = vec_to_point(x₀,point_guess)

        push!(branch.points, pd_corrected)
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

        if δ < δmin
            println("Stepsize to small: δ = $δ")
            break
        end

        # check if point is within bounds
        if parameterbounds !== nothing 
            if !all(>(0),(parameterbounds.min .< branch.points[end].parameters .< parameterbounds.max))
                break
            end
        end

    end
end


function continue!(branch)

    f = branch.f
    df = branch.df
    V = branch.tangents[end]
    δ=branch.δ
    MaxNumberofSteps=branch.MaxNumberofSteps
    δmin=branch.δmin
    δmax=branch.δmax
    parameterbounds = branch.parameterbounds
    con_par = branch.con_par
    NumberOfFails = branch.NumberOfFails

    fails = 0

    iterations=0
    for step = 1:MaxNumberofSteps

        # step != 1 || length(branch.points) != 0 ?  x₀ += δ₀*V : x₀
        point_guess = branch.points[end]
        # point_prev = branch.points[end]
        x₀ = vec(point_guess,con_par)
        # if step > 1
        #     point_prev = branch.points[end-1]
        #     # x₁ = vec(point_prev,con_par)
        #     # V = x₀ - x₁
        #     # @show V /= norm(V)
        # end
        x₀ += δ*V # new guess
        # @show display(δ)
        # @show display(V)
        x₀, iterations, V, converged  = newton(f, df, x₀, V, branch.points[end]; tol=1e-8)

        if ~converged
            println("No convergence")
            fails += 1
            if fails < NumberOfFails
                δ /= 2.0
                # branch = @set branch.points = branch.points[1:end-1] # remove last point
                continue
            else
                break
            end
        end
        fails = 0

        iterations < 5 ? δ *= 1.2 : δ /= 1.2
        δ = min(δ,δmax)

        point_corrected = vec_to_point(x₀,point_guess,con_par)

        push!(branch.points, point_corrected)
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

        if δ < δmin
            println("Stepsize to small: δ = $δ")
            break
        end

        # check if point is within bounds
        if parameterbounds !== nothing 
            if !all(>(0),(parameterbounds.min .< branch.points[end].parameters .< parameterbounds.max))
                println("Boundary hit.")
                break
            end
        end

    end
end

