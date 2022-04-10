## Processing

Since Dec 2, 2021, there is a single script "`full_process.py`" for processing leica (lif), nikon (nd2), and tiffs. It automatically corrects for motion, and phase (for bidirectional scans), and performs the usual sigmentation into ROIs. 

The script needs to be run in a terminal. To open the terminal from your [ctn landing page](https://ctn.physiologie.meduniwien.ac.at), click `New` and choose terminal.

Then, you can process your recording by issuing for example:  
`/data/useful/scripts/full_process.py "/data/example/full_path_to_your_file.nd2"`

For lif and nd2 files, the script will attempt to read important metadata, such as frequency from the file itself. If it can read also the pixel size, the script will automatically determine the size of spatial filters, which need to be approximately the diameter of cells (in pixels).
In case some of this does not work, or you know the metadata is not correct you can override the automatic choises by using optional arguments `--frequency` and/or `--spatial-filter`. 

For example, the following command will process the recording assuming it was done at 10Hz and will use three spatial filters 7,8, and 10:  
`/data/useful/scripts/full_process.py "/data/example/full_path_to_your_file.nd2" --frequency=10 --spatial-filter="7,8,10"`

For tifs, I did not implement automatic way of reading metadata, so you need to provide both of these arguments.


`--series` (for lifs):

Leica files often contain more than one recording ("series"). By default, the script assumes all the series belong together, and will try to stitch them.

If that does not work, or you wish to process only some series from the file, you can specify them using the argument `--series`. 

For example, `.../full_process.py "/full/path/to_file.lif" --series="Series011"`, will proess `Series011` from the file, and  `.../full_process.py "/full/path/to_file.lif" --series="Series011-13"` will first attempt to stitch `Series011`, `Series012` and `Series013`, and then process the resulting recording.

## Understanding output

In the default case, the script produces 6 output files: 
    
  - 2 videos in mp4 (original and corrected)
  - 2 plots in png (for phase shifts, and motion shifts)
  - 2 regions files for two different sizes of spatial filters (`.pkl` at the time of writing this tutotial, but we have intention to change to use other file format)
  
and saves them in a folder named for example `"/original/path/to_recording.lif_analysis/series_name/"`, where "`series_name`" is the name of the series for lif files (if you provided one), and "`all`" in all other cases.

## Other optional arguments

This document may not always reflect all the available options. You can always check the exact options and their succint description by running `/data/useful/scripts/full_process.py --help`. At the bottom of this document is the output of this command at the time of writing this text.

### Restricting to a concrete time window
  
If you want to restrict the analysis to only a part of the movie, from say 200th to 1234th second, you can use the `restrict` variable (`--restrict="200_1234"`).  
  
Restrict will also take the negative value for the end time point. For example, if you want to process all except the last 10s, you can enter `--restrict="0_-10"`.

When you use this option, the folder name inherits this argument as an additional suffix. For example, the folder name would be `"/original/path/to_recording.lif_analysis/series_name_0_-10s/"`

### I get too few ROIs!
  
By default, not all pixels of the recording will be taken to perform segmentation into ROIs, but only those which stick enough out of the background. A parameter that regulates this is called `--pixels-cutoff`. It's default value (0.02 at the time of writing this) is chosen so as to provide reasonable number of ROIs for the recordigs of mouse tissue slices. For other experiments it is possible that the value should be lower, or even put to 0. (If there is lot's of empty space, 0 might result in way too much ROIs, in which case, I suggest using a small positive number, e.g. 0.001.)

If you wish to try this option, add `--pixels-cutoff=0.001` to your command.

### Specifying a channel

By default, the processing is using the first channel (with index `0`). If you want to infer ROIs based on another channel, you can specify it by adding argument for example `--channel=1`, which will use the channel with index `1`. The output folder then get additinal suffix "`_c1`", e.g. `"/original/path/to_recording.lif_analysis/series_name_0_-10s_c1/"`

### More output

If you wish to follow more verbosely the processing, add `--verbose` to the command

### Notifications on slack

If you wish to get notifications on slack on the progress of processing, add `--notify` to the command. If this option does not work for you, we need to add your username manually to the script.

## Help from the terminal

`$full_process.py --help`

```
usage: full_process.py [-h] [--series SERIES] [--restrict RESTRICT]
                       [--verbose] [--spatial-filter SPATIAL_FILTER]
                       [--frequency FREQUENCY] [--channel CHANNEL]
                       [--pixels-cutoff PIXELS_CUTOFF] [--notify]
                       [--notification-user SLACK_USERLIST] [--debug]
                       --recording

Process a recording: correct phase, correct for motion, segment and save.

positional arguments:
  --recording           path to the recording, can be tif, lif, or nd2

optional arguments:
  -h, --help            show this help message and exit
  --series SERIES, -ser SERIES
                        name of the series (relevant only for lif files)
  --restrict RESTRICT   restrict analysis to the time interval (in seconds!),
                        e.g. "0_100" will only process first 100 seconds of
                        the movie, and "200_-10" will process the time region
                        from 200s until 10s before the end.
  --verbose             toggle verbose output on
  --spatial-filter SPATIAL_FILTER, -sp SPATIAL_FILTER
                        produce roi pickles with exactly these filter sizes,
                        e.g. -sp="5" or -sp="5+6" to produce simple rois with
                        indicated sizes, or sp="5,5+6" to produce both 5 and
                        5+6. Default (None) will guess four sizes based on
                        pxSize in the metadata if there are.
  --frequency FREQUENCY, -fr FREQUENCY
                        Frequency (frame rate) of the recording. Default
                        (None) will try to get it from the metadata.
  --channel CHANNEL     specify a channel to be used for construction of rois
                        (for multi-channel recording). Default is 0
  --pixels-cutoff PIXELS_CUTOFF
                        specifies the cutoff value of what will be considered
                        noise in the filtered image. The lower the value, the
                        fewer pixels will be considered noise and discarded,
                        and the result will have more ROIs. If you wish to
                        capture more potential cells, put this parameter to 0
                        (or even negative, but I did not test that). The
                        default is chosen as a reasonable compromise for most
                        applications in recordings of tissue slices.
  --notify              Triggers notification on slack when the script
                        starts/finishes.
  --notification-user SLACK_USERLIST, -nu SLACK_USERLIST
                        List of users to notify with the slack notifications
  --debug               toggle debug mode
```


```python

```
