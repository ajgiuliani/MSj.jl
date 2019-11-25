using Documenter, MSJ

makedocs(
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://ajgiuliani.github.io/MSJ.jl/stable/",
        assets = ["assets/favicon.ico"],
        analytics = "UA-132913317-2",
    ),    
    #source  = "src",
    #build   = "build",
    clean   = false,
    #doctest = true,
    modules = [MSJ],
    highlightsig = true,
    sitename="MSJ.jl",
    authors = "Alexandre Giuliani.",

    pages = [
        "Home"            => "index.md",
        "Tutorials"       => "tutorial.md",
        "Manual"          => "manual.md",
        "References"      => "reference.md",
        "Miscellaneous"   => "misc.md",
    ],
    #strict = true,
)

deploydocs(
    repo = "github.com/ajgiuliani/MSJ.jl.git",
    target = "build",
#    devbranch = "dev",
#    branch = "gh-pages",
#    devurl = "dev",
#    versions = ["stable" => "v^", "v#.#", devurl => "dev"]
)
