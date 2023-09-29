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
    QC_score_list_min = [] 
    QC_score_list_mean = [] 
    QC_score_all_bases = []
    for i in range(len(decoded)):
        intensity_array_int = decoded[i] 
        quality_list = []
        for j in range(len(intensity_array_int)):
            quality = (np.array(intensity_array_int[j]).max())/(np.array(intensity_array_int[j]).sum()) 
            quality_list.append(quality)
        quality_list =  [x if math.isnan(x) else x for x in quality_list]
        QC_score_min = np.array(quality_list).min() 
        QC_score_mean = np.array(quality_list).mean() 
        QC_score_list_min.append(QC_score_min)
        QC_score_list_mean.append(QC_score_mean)
        QC_score_all_bases.append(quality_list)
    df = decoded.to_features_dataframe()
    df['quality_minimum'] = QC_score_list_min
    df['quality_mean'] = QC_score_list_mean
    df['quality_all_bases'] = QC_score_all_bases
    return df

def process_experiment(exp_path, 
                        output, 
                        register = False, 
                        register_dapi = False,
                        masking_radius = 15, 
                        threshold = 0.002, 
                        sigma_vals = [1, 10, 30], # min, max and number
                        decode_mode = 'PRMC',
                        normalization_method = 'MH' # or MD
                ):
    

    # create output folder if not exists
    if not os.path.exists(output):
        os.makedirs(output)
    
    # load experiment file
    experiment = Experiment.from_json(exp_path)
    all_fovs = list(experiment.keys())

    # from output, find the FOVs not processed 
    csv_files = sorted(os.listdir(output))
    try:
        fovs_done = list(pd.DataFrame(csv_files)[0].str.split('.',expand = True)[0])
    except KeyError:
        print('no FOVS done')
        fovs_done = []
    
    # specify the files not done
    not_done = sorted(set(all_fovs).difference(set(fovs_done)))
    
    for i in not_done:
        print('decoding '+i)
        decoded = ISS_pipeline(experiment[i], experiment.codebook, register, register_dapi, masking_radius, threshold, sigma_vals, decode_mode, normalization_method)
        df = QC_score_calc(decoded)
        df.to_csv(output + i +'.csv')

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
