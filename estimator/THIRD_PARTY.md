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

To update it: download `happly.h` from upstream, replace the file and
update the commit above.
