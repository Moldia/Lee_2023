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
from starfish.core.spots.DecodeSpots.trace_builders import build_spot_traces_exact_match
import pandas as pd
from starfish.types import Axes, TraceBuildingStrategies
import warnings
warnings.filterwarnings('ignore')
from starfish.image import ApplyTransform, Filter, LearnTransform, Segment
from starfish.spots import FindSpots, DecodeSpots, AssignTargets
from starfish.types import Axes, FunctionSource
import pprint
from starfish.types import Features, Axes



def add_codebook(experiment_json_doc):
    experiment_json_doc['codebook'] = "codebook.json"
    return experiment_json_doc


def make_codebook_json(output_dir,codebook_csv):
    """ convert color code matrix in csv to json format"""
    codebook_array = []
    with open(codebook_csv, "r") as f:
        for line in f:
            line = line.rstrip('\n').split(',')
            codewords = []
            for r, colorcode in enumerate(line[1:]):
                codewords.append({Axes.ROUND.value: r,
                                  Axes.CH.value: int(colorcode)-1,
                                  Features.CODE_VALUE:1})
            codebook_array.append({Features.CODEWORD:codewords, Features.TARGET: line[0]})
    codebook = Codebook.from_code_array(codebook_array)
    codebook_json_filename = "codebook.json"
    codebook.to_json(os.path.join(output_dir, codebook_json_filename))

def get_tilepos(tilepos_xy_csv):
    tilexy = np.ndarray((0,2))
    with open(tilepos_xy_csv, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(',')
            tilexy = np.vstack([tilexy, [np.double(line[0]), np.double(line[1])]])
    return tilexy

def make_spacetx_format(path, codebook_csv,filenames = ['Base_1_stitched',
                                    'Base_2_stitched',
                                    'Base_3_stitched',
                                    'Base_4_stitched',
                                    'Base_5_stitched'],
                                    tile_dim=2000, 
                                    pixelscale =  0.1625, 
                                    channels = ["DAPI", "Cy3", "Cy5", "AF750", "AF488"],
                                    DO_decorators = ["AF750", "Cy5", "Cy3", "AF488"],
                                    folder_spacetx = 'SpaceTX_format', 
                                    nuclei_channel = 1):

    class ISSTile2D(FetchedTile):
        def __init__(self, file_path, fov):
            self.file_path = file_path
            self.fov = fov

        @property
        def shape(self) -> Tuple[int, ...]:
            return SHAPE

        @property
        def coordinates(self) -> Mapping[Union[str, Coordinates], Union[Number, Tuple[Number, Number]]]:
            return {
                Coordinates.X: (tilexy[self.fov, 0]*pixelscale, (tilexy[self.fov, 0] + tilesz)*pixelscale),
                Coordinates.Y: (tilexy[self.fov, 1]*pixelscale, (tilexy[self.fov, 1] + tilesz)*pixelscale),
                Coordinates.Z: (0.0, 0.0),
            }

        def tile_data(self) -> np.ndarray:
            return imread(self.file_path)


    class ISS2DPrimaryTileFetcher(TileFetcher):
        def __init__(self, path):
            self.path = path

        def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
            return ISSTile2D(os.path.join(self.path, "{}-{}/tile{}.tif".format(filenames[r], CHORDER[ch]+1, fov+1)), fov)


    class ISS2DAuxTileFetcher(TileFetcher):
        def __init__(self, path):
            self.path = path

        def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
            #return ISSTile2D(os.path.join(self.path, self.prefix + "{}.tif".format(fov+1)), fov)
            return ISSTile2D(os.path.join(self.path, "{}-{}/tile{}.tif".format(filenames[r], nuclei_channel ,fov+1)), fov)

    
    CHORDER = [channels.index(i) for i in DO_decorators]

    tilepos_xy_csv = path+'/preprocessing/ReslicedTiles/tilepos.csv'

    tilesz = tile_dim 
    SHAPE = {Axes.Y: tilesz, Axes.X: tilesz}
    num_tiles = len(os.listdir(path+'/preprocessing/ReslicedTiles/Base_1_stitched-1'))

    input_dir = path +'/preprocessing/ReslicedTiles'
    output_dir = path + '/' +folder_spacetx

    if not os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir)
        except:
            os.makedirs(output_dir)

    tilexy = get_tilepos(tilepos_xy_csv)


    write_experiment_json(
        path=output_dir, fov_count=num_tiles, tile_format=ImageFormat.TIFF,
        primary_image_dimensions={
            Axes.ROUND: len(filenames),
            Axes.CH: len(DO_decorators),
            Axes.ZPLANE: 1,
        },
        aux_name_to_dimensions={
            'nuclei': {
                Axes.ROUND: len(filenames),
                Axes.CH: 1,
                Axes.ZPLANE: 1,
            },
        },
        primary_tile_fetcher=ISS2DPrimaryTileFetcher(input_dir),
        aux_tile_fetcher={
            'nuclei': ISS2DAuxTileFetcher(input_dir),
        },
        postprocess_func=add_codebook,
        default_shape=SHAPE
    )

    make_codebook_json(output_dir,codebook_csv)

    from shutil import copyfile
    ls = os.listdir(output_dir)
    os.mkdir(os.path.join(output_dir, "originaljsons"))
    for file in ls:
        if file[-4:] == 'json':
            copyfile(os.path.join(output_dir, file), os.path.join(output_dir, "originaljsons", file))
            with open(os.path.join(output_dir, "originaljsons", file), 'r') as fr:
                with open(os.path.join(output_dir, file), 'w') as fw:
                    for line in fr:
                        fw.write("%s" % line.replace(output_dir.replace('\\', '\\\\') + '\\\\', ''))

