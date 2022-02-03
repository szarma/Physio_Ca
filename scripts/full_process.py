#!/opt/tljh/user/envs/physio_nc/bin/python
import os
import bioformats
import tifffile
import argparse
import islets
import pandas as pd
import numpy as np
import getpass
import sys
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Token for the bot, which will be used to post in slack for notifications
SLACK_BOT_TOKEN = 'xoxb-1186403589973-2442412187872-CtxhZTgETZGcfXGQ5uF2Y4kF'

# Member IDs of the users connected to slack
SLACK_USER_IDS = {
    'jupyter-johannes': 'U015WBY6A6M',
    'jupyter-sandra': 'U015MNRBCPM',
    'jupyter-marjan': 'U015J2J3NG2',
    'jupyter-nastja': 'U01C6HAJXNZ',
    'jupyter-srdjan': 'U01545E6T9V'
}


# noinspection PyBroadException
def create_movie_from_tif(recording, frequency, restrict):
    if frequency is None:
        raise ValueError("For tif files you need to specify the frequency.")
    else:
        try:
            fr = float(frequency)
        except ValueError:
            raise ValueError("Can not parse the frequency argument.")
    if restrict is None:
        slice_ = slice(None)
    else:
        try:
            t0, te = restrict.split("_")
            n0 = int(fr * float(t0))
            ne = int(fr * float(te))
            slice_ = slice(n0, ne)
        except ValueError:
            raise ValueError(f"Cannot parse the restrict argument")
    try:
        images = tifffile.memmap(recording, )[slice_]
    except:
        try:
            images = islets.utils.load_tif(recording, subindices = slice_)
        except:
            raise ImportError(f"Could not import the recording ({recording}).")

    movie = islets.cmovie(images, fr = fr)
    movie.filename = images.filename

    metadata = pd.Series({
        "Frequency": movie.fr,
        "SizeT": movie.shape[0],
        "SizeX": movie.shape[1],
        "SizeY": movie.shape[2],
        "Name": "all",
    })
    if restrict is not None:
        metadata["Name"] = metadata["Name"] + "_" + restrict

    try:
        metadata["frame_range"] = (n0, ne)
    except:
        pass

    return movie, metadata


def create_movie(recording, frequency, restrict, series, input_type, channel, verbose=True):
    try:
        if args.debug:
            global rec
    except:
        global rec
    if input_type == "tif":
        movie, metadata = create_movie_from_tif(recording, frequency, restrict)

    else:
        ##### input_type is either nikon or leica
        try:
            bioformats.javabridge.run_script("2+2")
        except:
            bioformats.javabridge.start_vm(class_path = bioformats.JARS)
            islets.utils._init_logger()

        rec = islets.Recording(recording)

        if input_type == "leica":
            if series is None or series=="all":
                from islets.Recording import parse_leica
                try:
                    parse_leica(rec, merge = True, verbose = verbose)
                    series = "all"
                except:
                    raise ValueError("""This appears to be a leica file, which can contain multiple series.
                    You did not specify the series, so I assumed you want to process all 
                    the series from the file, merged. I did not manage to do that. Try 
                    specifying the exact series, or the series range.""")
        else:
            if series is not None:
                raise UserWarning("For nikon files, the series argument is ignored.")
            series = "all"
        if restrict is None:
            restrict_ = None
        else:
            t0, te = restrict.split("_")
            t0 = float(t0)
            te = float(te)
            restrict_ = t0, te
        rec.import_series(series, restrict = restrict_, channel = channel, mode = "tif")
        metadata = rec.Series[series]["metadata"]
        images = rec.Series[series]["data"]

        if restrict is not None:
            metadata["Name"] = metadata["Name"] + "_" + restrict
        if channel > 0:
            metadata["Name"] = metadata["Name"] + "_c%i" % channel

        if frequency is None:
            fr = metadata["Frequency"]
        else:
            try:
                fr = float(frequency)
                metadata["Frequency"] = fr
            except:
                raise ValueError("Can not parse the frequency argument.")

        movie = islets.cmovie(images, fr = fr)
        movie.filename = images.filename
    del rec
    return movie, metadata


def slack_notify(text):
    client = WebClient(token = SLACK_BOT_TOKEN)
    try:
        user = []
        if args.slack_userlist is None:
            user = [getpass.getuser()]
        else:
            user = args.slack_userlist if type(args.slack_userlist) == list else [args.slack_userlist]
        for slack_user in user:
            if slack_user in SLACK_USER_IDS:
                client.chat_postMessage(channel = SLACK_USER_IDS[slack_user],
                                        text = text
                                        )
            else:
                print(f'User {slack_user} does not exist in slack dictionary. Message will not be sent.')
    except SlackApiError as e:
        assert e.response['error']
        print(f'Error while sending slack notification: {e.response["error"]}')

def main(args):
    if args.recording.lower().endswith("tif") or args.recording.lower().endswith("tiff"):
        input_type = "tif"
    elif args.recording.lower().endswith("lif"):
        input_type = "leica"
    elif args.recording.lower().endswith("nd2"):
        input_type = "nikon"
    else:
        raise NotImplementedError("Filetype not recognized. Currently, only tif, nd2, and lif are supported.")

    if args.notify:
        process = " ".join(sys.argv[1:])
        text = f"__Started__ the process: '{process}'"
        slack_notify(text = text)

    ## import data
    movie, metadata = create_movie(args.recording, args.frequency, args.restrict, args.series, input_type, args.channel,
                                   verbose = args.verbose)
    if movie[::len(movie)//10+1].std()/movie[::len(movie)//10+1].mean()<.1:
        raise NotImplementedError("Sorry, this movie seems to be just a still image.")
    ### prepare outputs
    outputDir = os.path.join(args.recording + "_analysis", metadata["Name"])
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)
        os.chown(outputDir, uid = -1, gid = 1002)
        os.chmod(outputDir, mode = 0o775)
    baseName = os.path.split(args.recording)[1] + "_" + metadata["Name"] + "_"
    baseName = os.path.join(outputDir, baseName)

    mp4Filename_original = baseName + "original.mp4"
    mp4Filename_corrected = baseName + "corrected.mp4"
    tifFilename_corrected = baseName + "corrected.tif"

    ####### Write the original movie #######
    ## because, for lif and nd2, it will be changed
    islets.utils.saveMovie(movie, filename = mp4Filename_original)

    ###### define or parse the spatial filters
    if args.spatial_filter is None:
        if input_type == 'tif':
            raise ValueError("For tif files, you need to specify the spatial filter. See help for more info.")
        else:
            try:
                filtSizes = islets.utils.get_filterSizes(metadata.pxSize)
                cell_half_width_in_px = int(np.ceil(6 / metadata.pxSize))
            except:
                raise ValueError(
                    "Could not read pixel size from the file's metadata. Please enter spatial filter(s) as an "
                    "argument. See help for more info")
    else:
        filtSizes = args.spatial_filter.split(",")
        filtSizes = [eval(el.replace("+", ",")) if "+" in el else (int(el),) for el in filtSizes]
        cell_half_width_in_px = min(sum(filtSizes, (np.inf,)))

    ####### Correcting ####################
    max_dev_rig = cell_half_width_in_px // 2 + 1
    if input_type == "tif":
        m_correct = islets.cmovie(
            tifffile.memmap(
                tifFilename_corrected,
                shape = movie.shape,
                dtype = "float32",
                photometric = "minisblack"
            ),
            fr = movie.fr
        )
    else:
        m_correct = movie

    # if args.debug: raise InterruptedError()

    # correct phase
    phase_shifts = islets.utils.correct_phase(movie, m_correct, max_dev = max_dev_rig, plot_name = baseName + "phase_shifts.png")
    np.savetxt(baseName+"phase_shifts.txt", phase_shifts, fmt = "%.4f", delimiter = ",")

    # correct for motion
    motion_shifts = islets.utils.motion_correct(movie, m_correct, max_dev = (max_dev_rig, max_dev_rig),
                                plot_name = baseName + "motion_shifts.png", verbose = args.verbose, mode = "full")
    np.savetxt(baseName+"motion_shifts.txt", motion_shifts, fmt = "%.4f", delimiter = ",")

    ####### Write corrected movie #######
    islets.utils.saveMovie(m_correct, filename = mp4Filename_corrected)

    ###### Start processing ######
    statistics = islets.utils.getStatImages(m_correct)

    for spFilt in filtSizes:
        if args.verbose: print("\t" * 2, "#" * 5, spFilt)
        pklBase = ".".join(map(str, spFilt))
        pickleFile = os.path.join(outputDir, pklBase + "_rois.pkl")
        if args.verbose: print("processing with filter size of ", spFilt)
        regions = islets.Regions(statistics, gSig_filt = spFilt, img_th = args.pixels_cutoff)
        if args.verbose:
            print(f"initialized with {len(regions.df)} rois.")
        regions.merge_closest(verbose = args.verbose)
        regions.sortInOrder()
        regions.metadata = metadata
        regions.calcTraces(m_correct)
        if "frame_range" in metadata and "Frequency" in metadata:
            regions.time += metadata.frame_range[0] / metadata.Frequency
        islets.utils.saveRois(regions, outputDir, filename = pklBase, add_date = False, formats = ["vienna"])
    ## Write the protocol template
    with open(baseName+"protocol.txt","w") as f:
        f.write("compound,concentration,begin,end\n,,,")

    ## Cleanup and notification
    bioformats.javabridge.kill_vm()

    if args.notify:
        text = text.replace("Started", "Finished")
        slack_notify(text)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'Process a recording: correct phase, correct for motion, segment and save.')
    parser.add_argument('recording', type = str,
                        help = 'path to the recording, can be tif, lif, or nd2')
    parser.add_argument('--series', '-ser', type = str, default = None,
                        help = 'name of the series (relevant only for lif files)')
    parser.add_argument('--restrict', type = str, default = None,
                        help = 'restrict analysis to the time interval (in seconds!), e.g. "0_100" will only process '
                               'first 100 seconds of the movie, and "200_-10" will process the time region from 200s '
                               'until 10s before the end.')
    parser.add_argument('--verbose', const = True, default = False, action = "store_const",
                        help = 'toggle verbose output on')
    # parser.add_argument('--rebin', default = None, type = "int", help = '''spatial rebin. By default, for lif and
    # nd2 files, the optimal rebin is inferred from the pxSize, if available. If it is not available, ''')
    parser.add_argument('--spatial-filter', "-sp", default = None,
                        help = '''produce roi pickles with exactly these filter sizes,
                        e.g. -sp="5" or -sp="5+6" to produce simple rois with indicated sizes,
                        or sp="5,5+6" to produce both 5 and 5+6. Default (None) will guess four
                        sizes based on pxSize in the metadata if there are.''')
    parser.add_argument('--frequency', "-fr", default = None, type = float,
                        help = '''Frequency (frame rate) of the recording. Default (None) will try to get it from the 
                        metadata.''')
    # parser.add_argument('--line-scan', default="none", type=str,
    #                     help='indicate if it is a line scan, and if yes, what kind ("single" or "multi")')
    parser.add_argument('--channel', default = 0, type = int,
                        help = 'specify a channel to be used for construction of rois (for multi-channel recording). '
                               'Default is 0')
    parser.add_argument('--pixels-cutoff', default = 0.02, type = float,
                        help =
                     """specifies the cutoff value of what will be considered noise in the filtered image. The lower
                        the value, the fewer pixels will be considered noise and discarded, and the result will have
                        more ROIs. If you wish to capture more potential cells, put this parameter to 0 (or even 
                        negative, but I did not test that). The default (0.02) is chosen as a reasonable compromise for
                        most applications in recordings of tissue slices."""
                        )
    parser.add_argument('--notify', dest = 'notify', action = 'store_true',
                        help = 'Triggers notification on slack when the script starts/finishes.')
    parser.add_argument('--notification-user', '-nu', default = None, dest = 'slack_userlist',
                        help = 'List of users to notify with the slack notifications')
    parser.add_argument('--debug', action = "store_true", dest = "debug",
                        help = 'toggle debug mode')

    args = parser.parse_args()

    if args.debug:
        args.verbose = True
        for k in args.__dict__.keys():
            print("%20s" % k, args.__dict__[k])

    try:
        main(args)
    except Exception as e:
        print (e)
    finally:
        sys.exit(0)
