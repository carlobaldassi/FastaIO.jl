# FastaIO.jl

| **Documentation**                                                         | **Build Status**                                             |
|:-------------------------------------------------------------------------:|:------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url][![][codecov-img]][codecov-url] |

Utilities to read/write FASTA format files in [Julia].

## Installation and usage

### Installation

To install the module, use Julia's package manager: start pkg mode by pressing <kbd>]</kbd> and then enter:

```
(v1.3) pkg> add FastaIO
```

Dependencies will be installed automatically.
The module can then be loaded like any other Julia module:

```
julia> using FastaIO
```

### Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**DEV**][docs-dev-url] &mdash; *in-development version of the documentation.*

See also the examples in the `examples/` directory.

[Julia]: http://julialang.org

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://carlobaldassi.github.io/FastaIO.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://carlobaldassi.github.io/FastaIO.jl/dev

[travis-img]: https://travis-ci.com/carlobaldassi/FastaIO.jl.svg?branch=master
[travis-url]: https://travis-ci.com/carlobaldassi/FastaIO.jl

[codecov-img]: https://codecov.io/gh/carlobaldassi/FastaIO.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/carlobaldassi/FastaIO.jl
