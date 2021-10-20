![Screenshot](docs/images/multipy-readme.svg)

------

## [Web documentation](https://multipy.readthedocs.io/en/latest/index.html)

------

## Installation

Follow these instructions:

```
git clone https://gitlab.multiscale.utah.edu/kamila/multipy.git
cd multipy
python setup.py install
```

If the installation was successful, you should be able to run the unit tests:

```
python -m unittest discover
```

and Python should now allow you to:

```
import multipy
```

------

## Local documentation build

To build the documentation locally first install the required libraries:

```
pip install sphinx
pip install jupyter-sphinx-theme
```

and then build the documentation:

```
cd docs
sphinx-build -b html . builddir
make html
```

To view the documentation in your web browser:

```
open _build/html/index.html
```
