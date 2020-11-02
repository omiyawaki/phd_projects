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
    "fig:malr-mon-lat"
    "fig:dalr-mon-lat"
    "fig:mpi-r1-temp"
    "fig:mpi-r1-decomp-nhmid"
    "fig:mpi-r1-decomp-shmid"
    "fig:mpi-r1-decomp-nhpole"
    "fig:mpi-r1-decomp-shpole")
   (LaTeX-add-bibliographies
    "references"))
 :latex)

