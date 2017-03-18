# FastaIO.jl

| **Documentation**                                                               | **PackageEvaluator**                                            | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:---------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.5-img]][pkg-0.5-url] [![][pkg-0.6-img]][pkg-0.6-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

Utilities to read/write FASTA format files in [Julia].

## Installation and usage

To install the module, use Julia's package manager:

```
julia> Pkg.add("FastaIO")
```

Dependencies will be installed automatically.
The module can then be loaded like any other Julia module:

```
julia> using FastaIO
```

### Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

See also the examples in the `examples/` directory.


[Julia]: http://julialang.org

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://carlobaldassi.github.io/FastaIO.jl/stable
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://carlobaldassi.github.io/FastaIO.jl/latest

[travis-img]: https://travis-ci.org/carlobaldassi/FastaIO.jl.svg?branch=master
[travis-url]: https://travis-ci.org/carlobaldassi/FastaIO.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/18ijkex153jkubw5/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/carlobaldassi/fastaio-jl/branch/master

[codecov-img]: https://codecov.io/gh/carlobaldassi/FastaIO.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/FastaIO.jl

[pkg-0.5-img]: http://pkg.julialang.org/badges/FastaIO_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=FastaIO
[pkg-0.6-img]: http://pkg.julialang.org/badges/FastaIO_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=FastaIO
