
using Revise
using Documenter, PeriodicNormalizationDDEs

ENV["INTERACTIVE"] = "true"
ENV["INTERACTIVE"] = "false"
makedocs(sitename="Periodic Normalization for Delay Differential Equations",
         format=Documenter.HTML(
         size_threshold_ignore=["NeuralMassModel.md","NeuralMassModelStatic.md,ActiveControlSystem.md"],
         example_size_threshold=9209111))

# makedocs(sitename="Periodic Normalization for Delay Differential Equations",
# format=Documenter.LaTeX())
