import pkg_resources
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning, )
    from .movies import movie as cmovie
    # from .movies import load as cload
from . import EventDistillery
from .examine import examine
from .CrossfilterApp import crossfilterApp
from .LineScan import LineScan
from .linescanner import plot_heatmap, plot_trace, examine
from .numeric import rebin
from .PicklePicker import PicklePicker
from .Regions import Regions, load_regions
from .Recording import Recording, saveMovie, parse_leica
from .utils import saveRois, get_filterSizes
from .fitting import fit_spikes

__version__ = pkg_resources.get_distribution('islets').version
