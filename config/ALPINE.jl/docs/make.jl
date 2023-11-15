using Documenter
using ALPINE

makedocs(
    sitename = "ALPINE",
    format = Documenter.HTML(),
    modules = [ALPINE]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
