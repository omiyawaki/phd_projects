(TeX-add-style-hook
 "20200806-002_w005"
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
    "sec:org75b15d6"
    "sec:org7fea45e"
    "fig:org865427b"
    "fig:orgdc967fb"
    "fig:orgc342aab"
    "fig:org91c292c"
    "fig:orgc52acad"
    "fig:org2c059cb"
    "sec:org6ccb884"
    "fig:org5bcf633"
    "fig:org99dcb89"
    "fig:orgcf86695"
    "sec:orga4e0468"
    "fig:org0dd4041"
    "fig:orgf832265"
    "sec:org67bc814")
   (LaTeX-add-bibliographies
    "../../../../../../mnt/c/Users/omiyawaki/Sync/papers/references"))
 :latex)

