To build html and latex from doc/ directory where conf.py is:
$ sphinx-build -b html . ${html_build_dir}
$ sphinx-build -b latex . ${latex_build_dir}

When building html, LIBRA_Software_Manual.pdf in the doc/ directory will be copied to ${html_build_dir}/_downloads/.

To build pdf from latex file:
$ pdflatex ${latex_build_dir}/LIBRA_Software_Manual.tex

Requirements:
Python, and sphinx. Additionally, pdflatex executable.
