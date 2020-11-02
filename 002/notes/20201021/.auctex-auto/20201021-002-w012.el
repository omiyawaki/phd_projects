(TeX-add-style-hook
 "20201021-002-w012"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("geometry" "margin=1in")))
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
    "hyperref"
    "geometry")
   (LaTeX-add-labels
    "sec:org4bb1a84"
    "sec:org19eabf8"
    "fig:comp-sh-0"
    "fig:comp-sh-16"
    "fig:comp-nh-0"
    "fig:interp-sh-16"
    "fig:diff-sh-16"
    "fig:interp-nh-0"
    "fig:diff-nh-0"
    "sec:orgdf69088"
    "fig:mse-nhmid-45"
    "fig:dmse-nhmid-45"
    "fig:dra-nhmid-45"
    "fig:dstf-lo-nhmid-45"
    "fig:dmse-l-nhmid-45"
    "fig:mse-nhmid-85"
    "sec:orgd23aaa9")
   (LaTeX-add-bibliographies
    "../../../../../../mnt/c/Users/omiyawaki/Sync/papers/references"))
 :latex)

