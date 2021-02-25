import pkg_resources
import warnings
import importlib

# noinspection PyUnresolvedReferences
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning, )
    warnings.filterwarnings("ignore", category=RuntimeWarning, )
    if importlib.util.find_spec("caiman") is not None:
        from caiman import load as cload
        from caiman import movie as cmovie
    else:
        try:
            from .movies import movie as cmovie
            print("Managed to import cmovie from .movies")
            from .movies import load as cload
            print("Managed to import cload from .load. Perhaps we can depart from full caiman?")
        except:
            print("But, not all that's necessary. Do not depart from caiman yet.")
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
