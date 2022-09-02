# TDA-for-Sea-Ice-Percolation
Code for 'Topological Data Analysis for Percolation in Arctic Melt Pond Evolution'

Dependencies:
Ripser: a TDA package which computes persistence diagrams - see https://ripser.scikit-tda.org/en/latest/
Persim: Package for plotting persistence diagrams - see https://persim.scikit-tda.org/en/latest/index.html

To install, copy the following into terminal if working with conda:

```
conda install -c conda-forge ripser persim
```

or if working with pip:

```
pip install cython ripser persim
```

Conda should automatically install cython; if using pip we have to specify this.


Code files:

my_functions.py : helper functions used throughout files
model_binary_to_fractal.py : takes directory containing folders {Run i} of binary images saved as .txt (1=pond, 0=ice), and calculates the signed Euclidean distance transform, persistence diagrams, persistence images calculated by 'binning' (see [1]) and plots images visualizing these as well as thumbnails to display if you wanted to visualize persistence images with PCA in an interactive data visualization package like Bokeh.
model_binary_to_fractal.py : takes directory containing folders {Run i} of binary images saved as .txt, and calculates PH fractal dimension for them.
model_binary_to_boundary_fractal.py : above, but computes PH fractal dimension of the boundary


[1]: Moon, C., Mitchell, S. A., Heath, J. E., & Andrew, M. (2019). Statistical inference over persistent homology predicts fluid flow in porous media. Water Resources Research, 55, 9592– 9603. https://doi.org/10.1029/2019WR025171
