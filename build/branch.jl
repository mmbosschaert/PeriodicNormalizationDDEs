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

function reverse_branch!(branch)
    reverse!(branch.points)
    reverse!(branch.tangents)
    reverse!(branch.stepsizes)
    
    for (i,p) in enumerate(branch.tangents)
        branch.tangents[i] *= -1 
    end
end
