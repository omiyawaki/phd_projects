(TeX-add-style-hook
 "updates_on_R16_and_ZX19"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("geometry" "margin=1in") ("parskip" "parfill")))
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
    "hyperref"
    "geometry"
    "parskip")
   (LaTeX-add-labels
    "sec:org0e938e6"
    "fig:orgb97d032"
    "fig:org4bb1748"
    "sec:orgab6736d"
    "sec:org5107872"
    "sec:orga14851e"
    "eq:org951f6d2"
    "fig:orge395c42"
    "sec:org976cbae"
    "eq:org2c27d6d"
    "fig:org5ed57dd"
    "fig:org646e7c8"
    "sec:orgb33d625"
    "sec:org311f93c"
    "eq:orge13f04e"
    "eq:orgc9e9a41"
    "eq:orge95430d"
    "eq:orgdee5432"
    "fig:org97fb8e3"
    "fig:orgb08e4b0"
    "sec:orgb01e4b2"
    "fig:org08e9671"
    "sec:org52a8193"
    "fig:org07325b1"
    "fig:org0e86b61"
    "fig:org506f3e2")
   (LaTeX-add-bibliographies
    "../../../../../../mnt/c/Users/omiyawaki/Sync/papers/references"))
 :latex)

