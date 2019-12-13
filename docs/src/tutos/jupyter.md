# Jupyter notebooks
[Jupyter](https://jupyter.org/) notebooks are web based documents that may contains both codes, figures and other textual elements (such as equations, links, ...).  Jupyter notebook may be installed easily using Julia:
```
julia> using Pkg
julia> Pkg.add("IJulia")
```
When IJulia is installed, then a notebook may be launch like this:
```
julia> using IJulia
julia> notebook()
```
The `notebook()` function should launch a web browser from which a new notebook  may be started. On each entry of the notebook code, Markdown or text may be inserted. Each line of code may be executed and will eventually return a result.
In the following, the tutorials are given in `Jupyter notebook` form and can be viewed using [nbviewer](http://nbviewer.jupyter.org/).

