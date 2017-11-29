# fcHiggs

Codes to calculate the cross sections and branching ratios of Higgs bosons studied in [arxiv.org/abs/arXiv:1711.08930](arXiv:1711.08930).

## Requirements

* C++ compiler supporting C++14 features ([Clang](http://clang.llvm.org/cxx_status.html) > 3.4, [GCC](https://gcc.gnu.org/projects/cxx-status.html) > 4.9),

* [LHAPDF 6](http://lhapdf.hepforge.org/) and the `NNPDF23_lo_as_0130_qed` PDF data. The latter can be installed by

```
lhapdf install NNPDF23_lo_as_0130_qed
```

## Usage

Each executable shows the input parameters. For instance, running

* `./bin/pph_neutral`

without any input will show

```
Usage: pph_neutral <m_H (GeV)> <tan(beta)> <cos(alpha-beta)> [output]
```

If `[output]` is not set, the output will be shown in `stdout`.

* `./bin/pph_neutral 400 1.0 0.05`

```
pph_neutral: p p --> H
pph_neutral: E_{CM} = 14.000000 TeV
pph_neutral: m_H = 400.000000 GeV
pph_neutral: tan(beta) = 1.000000, cos(alpha-beta) = 0.050000
pph_neutral: integrating for cross section ...
pph_neutral: ... done.
pph_neutral: total cross section = 0.039918 +- 0.000790 pb
```
