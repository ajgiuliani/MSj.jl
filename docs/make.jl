using Documenter, MSj

makedocs(
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        #prettyurls = !("local" in ARGS),
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://ajgiuliani.github.io/MSj.jl/stable/",
        assets = ["assets/favicon.ico"],
        analytics = "UA-132913317-2",
    ),    
    #source  = "src",
    #build   = "build",
    #clean   = false,
    #doctest = true,
    modules = [MSj],
    highlightsig = true,
    sitename="MSj.jl",
    authors = "Alexandre Giuliani",

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
    repo = "github.com/ajgiuliani/MSj.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "dev",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"]
)
