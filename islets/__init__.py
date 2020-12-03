import pkg_resources

from .examine import examine
from .CrossfilterApp import crossfilterApp
from .LineScan import LineScan
from .linescanner import plot_heatmap, plot_trace, examine
from .numeric import rebin
from .PicklePicker import PicklePicker
from .Regions1 import Regions, load_regions, getPeak2BoundaryDF, getGraph_of_ROIs_to_Merge, mergeBasedOnGraph
from .Recording1 import Recording, saveMovie, parse_leica
from .utils import saveRois, get_filterSizes

__version__ = pkg_resources.get_distribution('islets').version