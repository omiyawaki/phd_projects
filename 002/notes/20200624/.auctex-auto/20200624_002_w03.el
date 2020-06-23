(TeX-add-style-hook
 "20200624_002_w03"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
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
    "sec:orgbce2332"
    "sec:orgdc9e3cc"
    "sec:org58bd427"
    "sec:orgc9cefe0"
    "sec:org7795033"
    "sec:org2efab71"
    "sec:orge5788bd"
    "fig:org7a8841e"
    "fig:org5d1a2c0"
    "fig:org343fee3"
    "sec:org7f5e928"
    "fig:org652d7a0"
    "fig:org1b9c61d"
    "fig:orgf6fa329"
    "fig:orge273189"
    "sec:org382cad9")
   (LaTeX-add-bibliographies
    "../../../../../Sync/papers/references"))
 :latex)

