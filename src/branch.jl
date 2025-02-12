struct Branch
  points
  tangents
  stepsizes
  f
  df
  parameterbounds
  δ
  δmin
  δmax
  MaxNumberofSteps
  con_par
  NumberOfFails
  specialpoints
end

function Branch(; points, tangents, stepsizes, f, df, parameterbounds, δ, δmin, δmax, MaxNumberofSteps, con_par, NumberOfFails, specialpoints)
    Branch(points, tangents, stepsizes, f, df, parameterbounds, δ, δmin, δmax, MaxNumberofSteps, con_par, NumberOfFails, specialpoints)
end

function reverse_branch!(branch)
    reverse!(branch.points)
    reverse!(branch.tangents)
    reverse!(branch.stepsizes)
    
    for (i,p) in enumerate(branch.tangents)
        branch.tangents[i] *= -1 
    end
end

function find_corrected_double_hopf_points(hopf_branches, jet, τs; 
                                           max_iter=100, tol=1e-12)
    # Container for all corrected double Hopf points
    double_hopf_points = []

    for hopf_branch in hopf_branches
        # Locate indices of double Hopf points within this branch
        ind_double_hopf = locate_double_hopf(hopf_branch.points)

        # Correct the identified double Hopf points
        for p in hopf_branch.points[ind_double_hopf]
            p_corrected = LocateDoubleHopf(
                jet, p, τs;
                MaxIter = max_iter, 
                tol     = tol
            )
            push!(double_hopf_points, p_corrected)
        end
    end

    # Remove near-duplicate points
    PeriodicNormalizationDDEs.delete_one_from_approx_equal!(double_hopf_points)

    return double_hopf_points
end

"""
    find_genh_points(branch::Branch, jet, τs, par_indx; max_iter=100, tol=1e-10)

Identify and locate Generalized Hopf (genh) points for a single branch of type `Branch`.
"""
function find_genh_points(branch::Branch, jet, τs, par_indx; max_iter=100, tol=1e-10)
    # 1. Detect indices of GenH points in this branch
    ind_genh = detect_genh(jet, branch.points, τs)

    # 2. Refine the GenH points
    p_correct = [
        locate_genh(jet, branch, idx, par_indx, τs; MaxIter=max_iter, tol=tol) 
        for idx in ind_genh
    ]

    # 3. Filter out `nothing`
    p_correct = filter(!isnothing, p_correct)

    return p_correct
end

"""
    find_genh_points(branches::AbstractVector{Branch}, jet, τs, par_indx; max_iter=100, tol=1e-10)

Identify and locate Generalized Hopf (genh) points for multiple branches of type `Branch`.
"""
function find_genh_points(branches::AbstractVector{Branch}, jet, τs, par_indx; max_iter=100, tol=1e-10)
    # 1. Loop over branches, collecting GenH points from each
    all_points = Vector{Any}()

    for branch in branches
        branch_points = find_genh_points(branch, jet, τs, par_indx; max_iter=max_iter, tol=tol)
        push!(all_points, branch_points)
    end

    # 2. Flatten all points into a single vector
    return vcat(all_points...)
end

"""
    compute_nmfm_coefficients(double_hopf_points, jet, τs)

Compute normal form coefficients for each Double Hopf point by calling 
`normalform(jet, point, τs)`. If the computation fails for any point, 
a warning is printed and that particular point is skipped.

# Arguments
- `double_hopf_points`: A collection (e.g. `Vector`) of Double Hopf points.
- `jet`: The jet object or data required by `normalform`.
- `τs`: Additional parameters required by `normalform`.

# Returns
A new vector containing the computed normal form results for each point that
was successfully processed.
"""
function compute_nmfm_coefficients(double_hopf_points, jet, τs)
    new_points = DoubleHopf[]  # or specify a more precise type if known
    for point in double_hopf_points
        try
            push!(new_points, normalform(jet, point, τs))
        catch err
            @warn "Failed to calculate normal form coefficients for $(point)" exception=err
        end
    end
    return new_points
end
