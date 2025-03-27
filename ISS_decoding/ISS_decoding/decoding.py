import sys
import os
from typing import Mapping, Tuple, Union
import numpy as np
from starfish.types import Axes, Coordinates, Features, Number
from starfish import Codebook
from starfish.experiment.builder import FetchedTile, TileFetcher
from slicedimage import ImageFormat
from skimage.io import imread
from starfish.experiment.builder import write_experiment_json
#import matplotlib.pyplot as plt
import starfish
from copy import deepcopy
from starfish import data, FieldOfView, display, Experiment
from starfish.image import Filter
from starfish.spots import DetectPixels
from starfish.types import Axes, Features, Levels
from starfish import data, FieldOfView
from starfish.image import ApplyTransform, Filter, LearnTransform, Segment
from starfish.spots import FindSpots, DecodeSpots, AssignTargets
from starfish.types import Axes, FunctionSource, Levels
from starfish.core.expression_matrix.expression_matrix import ExpressionMatrix
from starfish.core.intensity_table.intensity_table import IntensityTable
test = os.getenv("TESTING") is not None
from starfish.core.spots.DecodeSpots.trace_builders import build_spot_traces_exact_match
import pandas as pd
from starfish.types import Axes, TraceBuildingStrategies
import warnings
warnings.filterwarnings('ignore')
from starfish.image import ApplyTransform, Filter, LearnTransform, Segment
from starfish.spots import FindSpots, DecodeSpots, AssignTargets
from starfish.types import Axes, FunctionSource
#import matplotlib
#import matplotlib.pyplot as plt
import pprint
from starfish.types import Features, Axes
#from starfish.util.plot import imshow_plane
test = os.getenv("TESTING") is not None
import math

def ISS_pipeline_dense(fov, codebook,
                 register = True, 
                 register_dapi = True,
                 masking_radius = 15,
                 threshold = 0.002,
                 sigma_vals = [1, 10, 30], # min, max and number
                 decode_mode = 'PRMC', # or MD
                 channel_normalization = 'MH' # if set to anything else, will do ClipPercentileToZero
                ):
    print('getting images')
    primary_image = fov.get_image(FieldOfView.PRIMARY_IMAGES) # primary images
    nuclei = fov.get_image('nuclei')
    dots = primary_image.reduce({Axes.CH, Axes.ROUND}, func="max") # reference round for image registration
    # register the raw image
    if register == True: 
        if register_dapi == True:
            nuclei_round1 = nuclei.sel({Axes.ROUND: 1, Axes.CH: 0, Axes.ZPLANE: 0})
            print('registering images based on nuclei stain')
            learn_translation = LearnTransform.Translation(reference_stack=nuclei_round1, axes=Axes.ROUND, upsampling=1000)
            transforms_list = learn_translation.run(nuclei)
        else:  
            print('creating reference images')
            dots = primary_image.reduce({Axes.CH, Axes.ROUND}, func="max") # reference round for image registration
            print('registering images')
            learn_translation = LearnTransform.Translation(reference_stack=dots, axes=Axes.ROUND, upsampling=100)
            transforms_list = learn_translation.run(primary_image.reduce({Axes.CH, Axes.ZPLANE}, func="max"))

        warp = ApplyTransform.Warp()
        registered = warp.run(primary_image, transforms_list=transforms_list, in_place=False, verbose=True)
        # filter registered data
        masking_radius = masking_radius
        filt = Filter.WhiteTophat(masking_radius, is_volume=False)
        #filtered = filt.run(primary_image, verbose=True, in_place=False)
        filtered = filt.run(registered, verbose=True, in_place=False)
    else: 
        print('not registering images')
        # filter raw data
        masking_radius = masking_radius
        filt = Filter.WhiteTophat(masking_radius, is_volume=False)
        filtered = filt.run(primary_image, verbose=True, in_place=False)

    # normalize the channel intensities
    print('normalizing channel intensities')
    if channel_normalization == 'MH':
        sbp = starfish.image.Filter.MatchHistograms({Axes.CH, Axes.ROUND})
    else: 
        sbp = starfish.image.Filter.ClipPercentileToZero(p_min=80, p_max=99.999, level_method=Levels.SCALE_BY_CHUNK)

    scaled = sbp.run(filtered, n_processes = 1, in_place=False)

    bd = FindSpots.BlobDetector(
        min_sigma=sigma_vals[0],
        max_sigma=sigma_vals[1],
        num_sigma=sigma_vals[2],
        threshold=threshold, # this is set quite low which means that we will capture a lot of signals
        measurement_type='mean',
    )

    # Instead of hardcoding channels, determine available channels dynamically.
    # Assuming primary_image is an xarray DataArray, we can get channel indices like this:
    channels = primary_image.coords[Axes.CH].values

    decoded_list = []
    intensities_list = []

    if decode_mode == 'PRMC':
        print('decoding with PerRoundMaxChannel')
        decoder = DecodeSpots.PerRoundMaxChannel(codebook=codebook)

    for ch in channels:
        # select a reference image for the current channel from round 0 and zplane 0
        channel_ref = primary_image.sel({Axes.ROUND: 0, Axes.CH: ch, Axes.ZPLANE: 0})
        print(f'locating spots for channel {ch}')
        spots = bd.run(reference_image=channel_ref, image_stack=scaled)

        if decode_mode == 'PRMC':
            decoded = decoder.run(spots=spots)
            print(f'Decoded spots: channel {ch}')
            # Calculate QC score for this channel's decoded spots and add to the list
            decoded_list.append(QC_score_calc(decoded))

        # Build spot traces for the current channel
        intensities = build_spot_traces_exact_match(spots)
        intensities_list.append(intensities)

    # Concatenate QC score dataframes for all channels if decoding was performed
    if decode_mode == 'PRMC' and decoded_list:
        decoded = pd.concat(decoded_list)
    else:
        decoded = None

    # Optionally, you might want to return or further process decoded, intensities_list, etc.
    return decoded

def ISS_pipeline(fov, codebook,
                register = True, 
                register_dapi = True,
                masking_radius = 15, 
                threshold = 0.002, 
                sigma_vals = [1, 10, 30], # min, max and number
                decode_mode = 'PRMC', # or MD
                channel_normalization = 'MH' # if set to anything else, will do ClipPercentileToZero
                ):

    print('getting images')
    primary_image = fov.get_image(FieldOfView.PRIMARY_IMAGES) # primary images
    nuclei = fov.get_image('nuclei')
    dots = primary_image.reduce({Axes.CH, Axes.ROUND}, func="max") # reference round for image registration
    # register the raw image
    if register == True: 
        if register_dapi == True:
            nuclei_round1 = nuclei.sel({Axes.ROUND: 0, Axes.CH: 0, Axes.ZPLANE: 0})
            print('registering images based on nuclei stain')
            learn_translation = LearnTransform.Translation(reference_stack=nuclei_round1, axes=Axes.ROUND, upsampling=1000)
            transforms_list = learn_translation.run(nuclei)
        else:  
            print('creating reference images')
            dots = primary_image.reduce({Axes.CH, Axes.ROUND}, func="max") # reference round for image registration
            print('registering images')
            learn_translation = LearnTransform.Translation(reference_stack=dots, axes=Axes.ROUND, upsampling=100)
            transforms_list = learn_translation.run(primary_image.reduce({Axes.CH, Axes.ZPLANE}, func="max"))
        
        warp = ApplyTransform.Warp()
        registered = warp.run(primary_image, transforms_list=transforms_list, in_place=False, verbose=True)
        # filter registered data
        masking_radius = masking_radius
        filt = Filter.WhiteTophat(masking_radius, is_volume=False)
        #filtered = filt.run(primary_image, verbose=True, in_place=False)
        filtered = filt.run(registered, verbose=True, in_place=False)
    else: 
        print('not registering images')
        # filter raw data
        masking_radius = masking_radius
        filt = Filter.WhiteTophat(masking_radius, is_volume=False)
        filtered = filt.run(primary_image, verbose=True, in_place=False)
        
    # normalize the channel intensities
    print('normalizing channel intensities')
    if channel_normalization == 'MH':
        sbp = starfish.image.Filter.MatchHistograms({Axes.CH, Axes.ROUND})
    else: 
        sbp = starfish.image.Filter.ClipPercentileToZero(p_min=80, p_max=99.999, level_method=Levels.SCALE_BY_CHUNK)

    scaled = sbp.run(filtered, n_processes = 1, in_place=False)
    
    bd = FindSpots.BlobDetector(
        min_sigma=sigma_vals[0],
        max_sigma=sigma_vals[1],
        num_sigma=sigma_vals[2],
        threshold=threshold, # this is set quite low which means that we will capture a lot of signals
        measurement_type='mean',
    )
    
    # detect spots using laplacian of gaussians approach
    dots_max = dots.reduce((Axes.ROUND, Axes.ZPLANE), func="max")
    print('locating spots')
    # locate spots in a reference image
    spots = bd.run(reference_image=dots_max, image_stack=scaled)
    
    if decode_mode == 'PRMC':
        print('decoding with PerRoundMaxChannel')
        decoder = DecodeSpots.PerRoundMaxChannel(codebook=codebook)
        decoded = decoder.run(spots=spots)
        # Build IntensityTable with same TraceBuilder as was used in MetricDistance

            
    elif decode_mode == 'MD':
        print('decoding with MetricDistance')
    # decode the pixel traces using the codebook
        decoder = DecodeSpots.MetricDistance(
            codebook=codebook,
            max_distance=1,
            min_intensity=1,
            metric='euclidean',
            norm_order=2,
            trace_building_strategy=TraceBuildingStrategies.EXACT_MATCH
        )
        decoded = decoder.run(spots=spots)
    
    intensities = build_spot_traces_exact_match(spots)

    # Get vector magnitudes, deal with empty tiles
    if intensities.size == 0: 
        print('No spots found')
    else:
        norm_intensities, vector_magnitude = codebook._normalize_features(intensities, norm_order=2)
    
    # Get distances
    distances = decoded.to_decoded_dataframe().data['distance'].to_numpy()

    return decoded

def QC_score_calc(decoded):
    # Lists to store quality metrics per image
    QC_score_list_min = [] 
    QC_score_list_mean = [] 
    QC_score_all_bases = []
    
    # Lists to store second peak ratio metrics per image
    second_peak_list_min = []
    second_peak_list_mean = []
    second_peak_all_bases = []
    
    # Helper function: compute second-peak ratio for one intensity array
    def compute_second_peak_ratio(intensity):
        arr = np.array(intensity)
        if arr.size == 0:
            return 0
        max_val = arr.max()
        if max_val == 0:
            return 0
        # Sort the array
        sorted_arr = np.sort(arr)
        # If only one element exists, return 0 for the second-peak
        if sorted_arr.size < 2:
            return 0
        second_val = sorted_arr[-2]
        return second_val / max_val
    
    # Loop over each image in decoded
    for i in range(len(decoded)):
        intensity_array_int = decoded[i] 
        quality_list = []
        second_peak_ratio_list = []
        
        # Loop over each cycle in the current image
        for j in range(len(intensity_array_int)):
            cycle = np.array(intensity_array_int[j])
            cycle_sum = cycle.sum()
            # Compute quality
            quality = cycle.max() / cycle_sum if cycle_sum != 0 else 0
            quality_list.append(quality)
            
            # Compute second-peak ratio for the cycle
            spr = compute_second_peak_ratio(intensity_array_int[j])
            second_peak_ratio_list.append(spr)
        
        # Remove any potential NaNs from the lists
        quality_list = [x if not math.isnan(x) else 0 for x in quality_list]
        second_peak_ratio_list = [x if not math.isnan(x) else 0 for x in second_peak_ratio_list]
        
        # Compute summary statistics for the current image
        QC_score_list_min.append(np.array(quality_list).min())
        QC_score_list_mean.append(np.array(quality_list).mean())
        QC_score_all_bases.append(quality_list)
        
        second_peak_list_min.append(np.array(second_peak_ratio_list).min())
        second_peak_list_mean.append(np.array(second_peak_ratio_list).mean())
        second_peak_all_bases.append(second_peak_ratio_list)
    
    # Assuming that decoded can produce a dataframe via to_features_dataframe()
    df = decoded.to_features_dataframe()
    df['quality_minimum'] = QC_score_list_min
    df['quality_mean'] = QC_score_list_mean
    df['quality_all_bases'] = QC_score_all_bases
    
    df['second_peak_ratio_min'] = second_peak_list_min
    df['second_peak_ratio_mean'] = second_peak_list_mean
    df['second_peak_ratio_all_bases'] = second_peak_all_bases
    
    return df


def process_experiment(exp_path, 
                        output, 
                        register=False, 
                        register_dapi=False,
                        masking_radius=15, 
                        threshold=0.002, 
                        sigma_vals=[1, 10, 30],  # min, max and number
                        decode_mode='PRMC',
                        normalization_method='MH',  # or other method
                        dense=False):
    # create output folder if not exists
    if not os.path.exists(output):
        os.makedirs(output)
    
    # load experiment file
    experiment = Experiment.from_json(exp_path)
    all_fovs = list(experiment.keys())

    # from output, find the FOVs not processed 
    csv_files = sorted(os.listdir(output))
    try:
        fovs_done = list(pd.DataFrame(csv_files)[0].str.split('.', expand=True)[0])
    except KeyError:
        print('no FOVS done')
        fovs_done = []
    
    # specify the files not done
    not_done = sorted(set(all_fovs).difference(set(fovs_done)))
    
    for i in not_done:
        print('decoding ' + i)
        if dense:
            # Override decode_mode if set to 'MD', because it is not applicable in dense mode.
            if decode_mode == 'MD':
                print('Warning: decode_mode "MD" is not applicable in dense mode. Overriding to "PRMC".')
                decode_mode_to_use = 'PRMC'
            else:
                decode_mode_to_use = decode_mode
            decoded = ISS_pipeline_dense(
                experiment[i], 
                experiment.codebook, 
                register=register, 
                register_dapi=register_dapi, 
                masking_radius=masking_radius, 
                threshold=threshold, 
                sigma_vals=sigma_vals, 
                decode_mode=decode_mode_to_use, 
                channel_normalization=normalization_method
            )
        else:
            decoded = ISS_pipeline(
                experiment[i], 
                experiment.codebook, 
                register=register, 
                register_dapi=register_dapi, 
                masking_radius=masking_radius, 
                threshold=threshold, 
                sigma_vals=sigma_vals, 
                decode_mode=decode_mode, 
                normalization_method=normalization_method
            )
        df = QC_score_calc(decoded)
        # Ensure proper file path concatenation
        output_path = os.path.join(output, i + '.csv')
        df.to_csv(output_path)

def concatenate_starfish_output(path, outpath,tag=''):
    
    import pandas as pd
    import os
    
    csv_files = sorted(os.listdir(path))
    print(len(csv_files))
    append_data = []
    for i in csv_files:
            print(i)
            df = pd.read_csv(path + i)
            df['fov'] = int(i.split('.')[0].split('_')[1])
            append_data.append(df)
    df_concat = pd.concat(append_data)
    spots_filt = df_concat[df_concat['target'].notna()]
    name=tag+'decoded.csv'
    df_concat.to_csv(os.path.join(outpath,name))
    return df_concat, spots_filt 

def plot_starfish_output(spots_file, 
                        dpi = 500, 
                        fig_size = (15,10), 
                        conversion = 0.1625, 
                        size_of_spots = 1):
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    mpl.rcParams['text.color'] = 'w'
    plt.style.use('dark_background')
    plt.rcParams["figure.figsize"] = fig_size
    mpl.rcParams['figure.dpi'] = dpi
    import pandas as pd
    
    df_concat = pd.read_csv(spots_file)
    spots_filt = df_concat[df_concat['target'].notna()]
    
    groups1 = spots_filt.groupby('target')
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    #io.imshow(dapi*10)
    for i, gene in enumerate(sorted(spots_filt.target.unique())):
        group1 = spots_filt[spots_filt.target == gene]
        ax.scatter(group1.xc/0.1625, group1.yc/0.1625, marker='.', linewidth=0, s=size_of_spots, label=gene)

    plt.gca().invert_yaxis()
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    import matplotlib.font_manager as fm
    fontprops = fm.FontProperties(size=10)
    scalebar = AnchoredSizeBar(ax.transData,
                               615.3846153846, '200 μm', 'lower right', 
                               pad=0.1,
                               color='white',
                               frameon=False,
                               size_vertical=5,
                               fontproperties=fontprops)
    ax.add_artist(scalebar)
    plt.axis('scaled')
    plt.axis('off')
    plt.title('Starfish decoding' + '\n' + 'Count: ' + str(spots_filt.shape[0]), size = 10)
    plt.show()
