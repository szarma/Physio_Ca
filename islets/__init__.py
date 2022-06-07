import warnings

import pkg_resources

# noinspection PyUnresolvedReferences
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning, )
    warnings.filterwarnings("ignore", category=RuntimeWarning, )
    # if importlib.util.find_spec("caiman") is not None:
    try:
        from caiman import movie as cmovie
    except ModuleNotFoundError:
        from .movies import movie as cmovie
        # print("Managed to import cmovie from .movies")
    # try:
    #     from caiman import load as cload
    # except:
    #     try:
    #         from .movies import load as cload
    #         print("Managed to import cload from .load. Perhaps we can depart from full caiman?")
    #     except:
    #         print("could not import load from anywhere.")

        # except:
        #     print("But, not all that's necessary. Do not depart from caiman yet.")
    from . import EventDistillery
    from .examine import examine
    from .CrossfilterApp import crossfilterApp
    from .LineScan import LineScan
    from .linescanner import plot_heatmap, plot_trace, examine
    from .numeric import rebin
    from .PicklePicker import PicklePicker
    from .Regions import Regions, load_regions
    from .Recording import Recording, parse_leica
    from .utils import saveRois, get_filterSizes, saveMovie
    from .fitting import fit_spikes
    from .protocol import Protocol

__version__ = pkg_resources.get_distribution('islets').version
