using Documenter, msJ

makedocs(
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://ajgiuliani.github.io/msJ.jl/stable/",
        assets = ["assets/favicon.ico"],
        analytics = "UA-132913317-2",
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
    #strict = true,
)

deploydocs(
    repo = "github.com/ajgiuliani/msJ.jl.git",
    target = "build",
    devbranch = "dev",
    branch = "gh-pages",
    devurl = "dev",
#    versions = ["stable" => "v^", "v#.#", devurl => devurl]
)
