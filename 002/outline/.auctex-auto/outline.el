(TeX-add-style-hook
 "outline"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "ametsocV5"
    "ametsocV510"
    "amsmath"
    "amsfonts"
    "amssymb"
    "bm"
    "mathptmx"
    "newtxtext"
    "newtxmath")
   (LaTeX-add-labels
    "eq:mse"
    "sec:results"
    "subsec:asym"
    "eq:mse-toasfc"
    "eq:mse-toasfc-approx"
    "eq:sfc"
    "eq:thf-full"
    "eq:thf-deep"
    "eq:d-deep"
    "eq:thf-shallow"
    "eq:in-phase"
    "eq:d-shallow"
    "fig:mpi-r1-ann"
    "fig:mpi-r1-dev"
    "fig:mpi-r1-decomp-mid"
    "fig:mpi-r1-decomp-pole"
    "fig:echam-rce"
    "fig:echam-rae"
    "fig:malr-mon-lat"
    "fig:dalr-mon-lat")
   (LaTeX-add-bibliographies
    "references"))
 :latex)
