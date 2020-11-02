(TeX-add-style-hook
 "20200924-002-w008"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "inputenc"
    "fontenc"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "amssymb"
    "capt-of"
    "hyperref")
   (LaTeX-add-labels
    "sec:orgc88dfe1"
    "sec:orgf4a5237"
    "fig:org3f41b82"
    "fig:org77bd9bf"
    "sec:orgb9af566"
    "fig:orgaf7c689"
    "fig:org2f02c30"
    "fig:org87d79ec"
    "sec:org792a68a"
    "fig:org734dc14"
    "sec:orgddffd35"
    "fig:org8bf6a4b"
    "fig:org7fca9bb"
    "sec:org4206401"
    "fig:orgfd6ae58"
    "sec:orgb3338f0"
    "fig:orgae48b0a"
    "sec:org70614aa"
    "fig:org1573ac9"
    "fig:org6023f92"
    "sec:orgb97cfc7"
    "fig:orgaa39c6c"
    "fig:orgf470bb8"
    "sec:orgad99e4c"
    "sec:org1d395c5"
    "sec:org50c378b"
    "sec:orgb162501")
   (LaTeX-add-bibliographies
    "../../../../../../mnt/c/Users/omiyawaki/Sync/papers/references"))
 :latex)

