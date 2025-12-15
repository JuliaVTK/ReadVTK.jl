using Documenter
import Pkg
using VTKBase
using ReadVTK

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(VTKBase, :DocTestSetup, :(using VTKBase); recursive = true)
DocMeta.setdocmeta!(ReadVTK, :DocTestSetup, :(using ReadVTK); recursive = true)

# Path to markdown file containing docs for VTKBase.jl
# We copy the file to src/external/
vtkbase_docs_src = joinpath(dirname(dirname(pathof(VTKBase))),  # VTKBase directory
                            "docs",
                            "src",
                            "VTKBase.md")
isfile(vtkbase_docs_src) || error("file not found: $vtkbase_docs_src")
vtkbase_docs = joinpath("external", "VTKBase.md")
cp(vtkbase_docs_src, joinpath(@__DIR__, "src", vtkbase_docs); force = true)

# Make documentation
makedocs(
         # Specify modules for which docstrings should be shown
         modules = [VTKBase, ReadVTK],
         # Set sitename to ReadVTK
         sitename = "ReadVTK.jl",
         # Provide additional formatting options
         format = Documenter.HTML(
                                  # Disable pretty URLs during manual testing
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  # Explicitly add favicon as asset
                                  # assets = ["assets/favicon.ico"],
                                  # Set canonical URL to GitHub pages URL
                                  canonical = "https://juliavtk.github.io/ReadVTK.jl/stable"),
         # Explicitly specify documentation structure
         pages = [
           "Home" => "index.md",
           "Reference" => "reference.md",
           "VTKBase.jl" => vtkbase_docs,
           "Contributing" => "contributing.md",
           "License" => "license.md",
         ])

deploydocs(repo = "github.com/JuliaVTK/ReadVTK.jl",
           devbranch = "main",
           # Only push previews if all the relevant environment variables are non-empty.
           push_preview = all(!isempty,
                              (get(ENV, "GITHUB_TOKEN", ""),
                               get(ENV, "DOCUMENTER_KEY", ""))))
