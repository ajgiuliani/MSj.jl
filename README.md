# msJ.jl

*A mass spectrometry package for Julia*

[![Build Status](https://travis-ci.org/ajgiuliani/msJ.jl.svg?branch=master)](https://travis-ci.org/ajgiuliani/msJ.jl)

[![Coverage Status](https://coveralls.io/repos/github/ajgiuliani/msJ.jl/badge.svg?branch=master)](https://coveralls.io/github/ajgiuliani/msJ.jl?branch=master)
[![codecov](https://codecov.io/gh/ajgiuliani/msJ.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ajgiuliani/msJ.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ajgiuliani.github.io/msJ.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ajgiuliani.github.io/msJ.jl/dev/)


## Installation
This package is unregistered. It can be installed either with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:
```julia
(v1.2) pkg> add https://github.com/ajgiuliani/msJ.jl
```
or using the package API:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/ajgiuliani/msJ.jl"))
```

## Documentation
Documentation is available [here](https://ajgiuliani.github.io/msJ.jl/dev).


## Project Status
The package is developed for Julia 1.0 and above on Linux and OSX. It should be working on Windows as well.


## Other packages
* [mzXML](https://github.com/timholy/mzXML.jl): Load mass spectrometry mzXML files.
