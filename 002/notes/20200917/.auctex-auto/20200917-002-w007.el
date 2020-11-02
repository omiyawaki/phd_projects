(TeX-add-style-hook
 "20200917-002-w007"
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
    "sec:org830e49f"
    "sec:org763d95b"
    "sec:org64e27bf"
    "fig:orgd386116"
    "fig:org196191e"
    "fig:org1360b06"
    "sec:orgb2afb46"
    "fig:org41c24b1"
    "sec:orgdd7ad82"
    "fig:org3c66fec"
    "fig:org0e1059c"
    "fig:org28f944a"
    "fig:org8e028a5"
    "fig:org660a9dd"
    "fig:org48b3b97"
    "sec:orgdc5106c")
   (LaTeX-add-bibliographies
    "../../../../../../mnt/c/Users/omiyawaki/Sync/papers/references"))
 :latex)

