# remove this function
function stability(M, points, τs; neigs=13, maxit=100)
  for (i, p) in enumerate(points)
    dep = DEP(M(repeat(p.coords, 1, length(τs)), p.parameters), vec(τs))
    try
      λ, _ = iar_chebyshev(dep, maxit=maxit, neigs=neigs)
      points[i] = @set p.stability = λ
    catch e
      println("An error occurred in iar_chebyshev of point number $i: ", e)
      points[i] = @set p.stability = points[i-1].stability
      continue
      # Handle the error or rethrow it
      # For instance, you can set a default value or leave the point unchanged
    end
  end
end

function stability(jet::NamedTuple, p, τs; neigs=13, maxit=100, σ=nothing)
  dep = DEP(jet.Ms(repeat(p.coords, 1, length(τs)), p.parameters), vec(τs))
  if σ === nothing
    λ, _ = iar_chebyshev(dep, maxit=maxit, neigs=neigs)
  else
    println("Using σ")
    λ, _ = iar_chebyshev(dep, maxit=maxit, neigs=neigs, σ=σ)
  end

  @set p.stability = λ
end

function stability(jet::NamedTuple, branch::NamedTuple, τs; neigs=13, maxit=100)
  for (i, p) in enumerate(branch.points)
    branch.points[i] = stability(jet, p, τs; neigs=neigs, maxit=maxit)
  end
end
