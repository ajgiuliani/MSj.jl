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
        "Home"                             => "index.md",
        "Tutorials"                        => Any[
            "The Julia language"           => "tutos/julia.md",
            "Jupyer notebooks"             => "tutos/jupyter.md",
            "MSj"                          => "tutos/MSj.md",
            ],
        "Manual"                           => Any[
            "Introduction"                 => "man/introduction.md",
            "Public elements"              => "man/public.md",
            "Data types"                   => "man/types.md",
            "File Information"             => "man/information.md",
            "Importing data"               => "man/importing.md",
            "Exporting data"               => "man/exporting.md",
            "Combining and filtering data" => "man/filtering.md",
            "Processing"                   => "man/processing.md",
            "Properties calculations"      => "man/calculations.md",
            "Plotting"                     => "man/plotting.md",
            ],
        "References"                       => "reference.md",
        "Miscellaneous"                    => "misc.md",
    ],
    #strict = true,
)

deploydocs(
    repo = "github.com/ajgiuliani/MSj.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    #devurl = "stable",
    versions = ["stable" => "v^", "v#.#"]
)
