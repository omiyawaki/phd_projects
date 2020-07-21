(TeX-add-style-hook
 "20200801_002_w005"
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
    "sec:orgb5ececc"
    "sec:orgaf98902"
    "sec:org63e75b4"
    "fig:org076a73b"
    "fig:org719411c"
    "sec:org1637c9e"
    "sec:org4370cc0"
    "sec:org55f55ff"
    "sec:orga0cd988"
    "sec:org95eef24"
    "sec:orgb10be3e"
    "fig:org5cdef74"
    "fig:org864fef8"
    "sec:org1329435"
    "tab:org64baad4"
    "tab:org4d3e2b8"
    "tab:org6fa0a98"
    "tab:org8748e47"
    "sec:org1d53cf8"
    "fig:org213f7f6"
    "sec:orgc56ad75"
    "fig:org84940c7"
    "fig:orgf87d332"
    "sec:org3161d01"
    "fig:org29f0494"
    "fig:org26b6917"
    "fig:orge633f0b"
    "sec:orgb04731d"
    "sec:org9386ade"
    "fig:org39dc232"
    "sec:orgfd7cab0"
    "fig:orgf2b3688"
    "sec:org358b4c0"
    "fig:org1cfa519"
    "sec:org9b1c8ac"
    "fig:orgea68b49"
    "sec:org1c6f734")
   (LaTeX-add-bibliographies
    "../../../../Sync/papers/references"))
 :latex)

