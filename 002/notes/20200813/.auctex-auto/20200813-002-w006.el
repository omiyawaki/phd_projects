(TeX-add-style-hook
 "20200813-002-w006"
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
    "sec:orgf4b5fd4"
    "sec:org202f2f8"
    "fig:org6bb4767"
    "sec:org54f98fe"
    "sec:orgf1e72ac"
    "sec:orgc858e38"
    "fig:org43fc374"
    "fig:org1fa1e68"
    "sec:org5643330"
    "sec:org61de612"
    "sec:org4e21439"
    "fig:orga08856b"
    "fig:org17106c3"
    "sec:orgece31ed"
    "sec:orgbbb8d7d"
    "sec:org817b123"
    "fig:org845a067"
    "fig:orgdb46e7a"
    "fig:org5617d48"
    "sec:org8bd0b2a"
    "sec:org5fe1813"
    "sec:org0a06548"
    "fig:org7193867"
    "fig:org075c5db"
    "fig:org80bb966"
    "sec:orgf9e0a67"
    "fig:org323adda"
    "sec:org172d899"
    "fig:org41a190c"
    "fig:org904dc58"
    "fig:org8cb1bab"
    "fig:orgf9db9e1"
    "fig:org0e7e926"
    "fig:org18c58bc"
    "fig:orga993cae"
    "sec:org65d6f89")
   (LaTeX-add-bibliographies
    "../../../../../../mnt/c/Users/omiyawaki/Sync/papers/references"))
 :latex)

