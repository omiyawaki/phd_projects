(TeX-add-style-hook
 "20200617_002_w002"
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
    "sec:orgfe30580"
    "sec:orge7a59bf"
    "sec:org9002395"
    "sec:org3eea186"
    "sec:org8d810ee"
    "sec:orgeb1a078"
    "sec:org8cf3576"
    "sec:orgcb2ee01"
    "sec:org3a05cf2"
    "fig:orga68063a"
    "fig:org04cdb55"
    "fig:org0e8e604"
    "sec:org4a486a4"
    "sec:orgc196c7a")
   (LaTeX-add-bibliographies
    "../../../../../Sync/Papers/references"))
 :latex)

