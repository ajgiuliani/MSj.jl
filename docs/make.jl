using Documenter, msJ

makedocs(
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        #canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
        #assets = ["assets/favicon.ico"],
        #analytics = "",
    ),
    #source  = "src",
    #build   = "build",
    clean   = false,
    #doctest = true,
    modules = [msJ],
    highlightsig = true,
    sitename="msJ.jl",
    authors = "Alexandre Giuliani.",

    pages = [
        "Home"           => "index.md",
        "Manual"         => "manual.md",
        "Tutorial"       => "tutorial.md",
        "Reference"      => "reference.md",
        "Miscellaneous"  => "misc.md",
    ],
    strict = true,
)

deploydocs(
    repo = "ithub.com/ajgiuliani/msJ.jl.git",
    target = "build",
)
