using Documenter
import Pkg
using ReadVTK

# Define module-wide setups such that the respective modules are available in doctests
DocMeta.setdocmeta!(ReadVTK, :DocTestSetup, :(using ReadVTK); recursive=true)

# Make documentation
makedocs(
    # Specify modules for which docstrings should be shown
    modules = [ReadVTK],
    # Set sitename to ReadVTK
    sitename="ReadVTK.jl",
    # Provide additional formatting options
    format = Documenter.HTML(
        # Disable pretty URLs during manual testing
        prettyurls = get(ENV, "CI", nothing) == "true",
        # Explicitly add favicon as asset
        # assets = ["assets/favicon.ico"],
        # Set canonical URL to GitHub pages URL
        canonical = "https://trixi-framework.github.io/ReadVTK.jl/stable"
    ),
    # Explicitly specify documentation structure
    pages = [
        "Home" => "index.md",
        "Reference" => "reference.md",
        "Contributing" => "contributing.md",
        "License" => "license.md"
    ],
    strict = true # to make the GitHub action fail when doctests fail, see https://github.com/neuropsychology/Psycho.jl/issues/34
)

deploydocs(
    repo = "github.com/trixi-framework/ReadVTK.jl",
    devbranch = "main",
    push_preview = true
)
