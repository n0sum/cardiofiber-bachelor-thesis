# Third-party dependencies

## happly

Single-header PLY reader/writer, by Nicholas Sharp.

- upstream : https://github.com/nmwsharp/happly
- license  : MIT
- commit   : 8a606309daaa680eee495c8279feb0b704148f4a

It was a git submodule in the previous repository. It has been **vendored**
(copied into `include/`) so that the repository builds without network access and without
depending on an upstream that can move: a necessary condition for
the build to be reproducible over time.

Local change: in the ASCII reader (`parseNext` of scalar and list properties) floating-point
tokens are parsed with `std::strtod`/`strtof` instead of `>>`, which with libstdc++ does not read
`nan`, `-nan`, `inf` and silently gives 0 (the non-estimable faces of format D are NaN); integers
are read with `>>` as upstream.

Local change: `#include <cstdint>` added, for `uint8_t` etc.; GCC 16 no longer includes it through
the other standard headers.

To update it: download `happly.h` from upstream, replace the file,
re-apply the local changes above and update the commit.
