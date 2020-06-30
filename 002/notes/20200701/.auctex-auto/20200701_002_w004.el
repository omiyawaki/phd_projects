(TeX-add-style-hook
 "20200701_002_w004"
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
    "sec:org223d7eb"
    "sec:org3a59a29"
    "sec:orga5fc947"
    "fig:org1588f96"
    "fig:org4a5a954"
    "sec:org8053936"
    "sec:org2467533"
    "sec:orgbd89dd0"
    "sec:orgb57f76e"
    "sec:orgdd37950"
    "sec:orga0d8db3"
    "fig:orgc1199f7"
    "fig:org30a1620"
    "sec:org819881b"
    "tab:orgf86a5a0"
    "tab:org04774cc"
    "tab:org31af5d2"
    "tab:orgaaee4f9"
    "sec:org63a1848"
    "fig:org40e5272"
    "sec:orgb2fa02d"
    "fig:org90bf693"
    "fig:org7296393"
    "sec:orga307dc7"
    "fig:orge3fdd15"
    "fig:org8e4c2fb"
    "fig:org6160dfc"
    "sec:org2f978e1"
    "sec:org32e6e8d"
    "fig:orgf8082a4"
    "sec:org3906601"
    "fig:org11c0aa0"
    "sec:orge1b7cb9"
    "fig:orgb8d6b32"
    "sec:org030ace5")
   (LaTeX-add-bibliographies
    "../../../../../Sync/papers/references"))
 :latex)

