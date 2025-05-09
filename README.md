# Grape extension to CepGen

A CepGen plugin for the [Grape](https://arxiv.org/abs/hep-ph/0012029v2) ep &rarr; e(&gamma;&gamma;&rarr;&#8467;<sup>+</sup>&#8467;<sup>-</sup>)p process definition.

## Installation and usage

A version of [CepGen](cepgen/cepgen) is required for all stages of the building and running of this extension.
Once set, the `CEPGEN_PATH` environment variable should be set to the installation path (possibly the location of your `cepgen` directory if directly cloned from the main repository, or the standard environment libraries path if installed from a RPM/DEB package.
Using the recommended building recipe of CepGen, this can be done the following way:

```bash
export CEPGEN_PATH=/path/to/your/cepgen/install
export LD_LIBRARY_PATH=$CEPGEN_PATH/lib64:$LD_LIBRARY_PATH
```

The Grape add-on shared library can then be built the usual "CMake way":

```bash
mkdir build && cd build
cmake [-GNinja] ..
[ninja|make -j]
```

This will provide you with a `libCepGenGrape.so` shared object than can be passed directly to the main CepGen executable, and steered using one of the example cards provided in the [cards](https://github.com/cepgen/CepGenGrapeProcess/tree/main/cards) directory, for instance:

```bash
$CEPGEN_PATH/bin/cepgen -a [path/to/]libCepGenGrape.so -i [path/to/cards/]grape_ep_epee_cfg.py
```
