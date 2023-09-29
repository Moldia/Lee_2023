from xml.dom import minidom
from tqdm import tqdm
import os
import pandas as pd
import tifffile
import numpy as np
import re
import shutil
from os.path import join
import ISS_processing.preprocessing as preprocessing
import ashlar.scripts.ashlar as ashlar
import cv2
import math
import mat73
import pathlib
import xml.etree.ElementTree as ET
from aicspylibczi import CziFile







def customcopy(src, dst):
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    shutil.copyfile(src, dst)
def zen_OME_tiff(exported_directory, output_directory, channel_split=2, cycle_split=1, num_channels=5):
    '''
    This function makes OME-TIFF files from files exported from as tiff from ZEN, through the process_czi or to the deconvolve_czi functions.
    
    Note: This function assumes that you are using the Nilsson SOP for naming files. It will work on 1-tile sections.
    Args:
    - exported_directory: directory containing exported TIFF files.
    - output_directory: directory to save the processed files.
    - channel_split, cycle_split: indices for splitting filenames.
    - num_channels: number of channels in the images.
    
    Returns:
    - None. Writes processed images to output_directory.
    '''

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Filter out TIFF files from the directory
    all_files = os.listdir(exported_directory)
    tiff_files = [file for file in all_files if '.tif' in file]

    # Split the filenames to extract tiles, channels, and rounds
    split_tiles_df = pd.DataFrame(tiff_files)[0].str.split('m', expand=True)
    split_channels_df = split_tiles_df[0].str.split('_', expand=True)
    tiles = list(np.unique(split_tiles_df[1]))
    channels = list(np.unique(split_channels_df[channel_split]))
    rounds = list(np.unique(split_channels_df[cycle_split]))

    # Iterate through rounds to process files
    for _, round_number in enumerate(rounds):
        tiff_files_round = [file for file in tiff_files if f'Base_{round_number}_' in file]
        metadata_files = [file for file in all_files if 'info.xml' in file]
        metadata_files_round = [file for file in metadata_files if f'_{round_number}_' in file]

        # Parse metadata XML files to extract tile positions
        for metadata_file in metadata_files_round:
            xml_doc = minidom.parse(os.path.join(exported_directory, metadata_file))
            tiles_xml, x_coords, y_coords = [], [], []
            bounds_elements = xml_doc.getElementsByTagName('Bounds')
            for elem in bounds_elements:
                tiles_xml.append(int(elem.attributes['StartM'].value))
                x_coords.append(float(elem.attributes['StartX'].value))
                y_coords.append(float(elem.attributes['StartY'].value))
                
            unique_tiles_xml = list(np.unique(tiles_xml))
            position_df = pd.DataFrame({
                'x': x_coords[:len(unique_tiles_xml)],
                'y': y_coords[:len(unique_tiles_xml)]
            })
            positions = np.array(position_df).astype(int)

        # Write processed images to OME-TIFF format
        with tifffile.TiffWriter(os.path.join(output_directory, f'cycle_{round_number}.ome.tif'), bigtiff=True) as tif:
            for i in tqdm(range(len(tiles))):
                position = positions[i]
                tile = tiles[i]
                tiff_files_tile = [file for file in tiff_files_round if f'm{tile}' in file and '._' not in file]
                stacked_images = np.empty((num_channels, 2048, 2048))

                for idx, image_file in enumerate(sorted(tiff_files_tile)):
                    image_data = tifffile.imread(os.path.join(exported_directory, image_file))
                    stacked_images[idx] = image_data.astype('uint16')

                pixel_size = 0.1625
                metadata = {
                    'Pixels': {
                        'PhysicalSizeX': pixel_size,
                        'PhysicalSizeXUnit': 'µm',
                        'PhysicalSizeY': pixel_size,
                        'PhysicalSizeYUnit': 'µm'
                    },
                    'Plane': {
                        'PositionX': [position[0] * pixel_size] * stacked_images.shape[0],
                        'PositionY': [position[1] * pixel_size] * stacked_images.shape[0]
                    }
                }
                tif.write(stacked_images.astype('uint16'), metadata=metadata)


def leica_mipping(input_dirs, output_dir_prefix, image_dimension=[2048, 2048]):
    """
    Process and MIP (maximum intensity projection) microscopy image files exported from Leica as TIFFs.

    Parameters:
    - input_dirs: List of file paths to the input directories.
    - output_dir_prefix: Prefix for the output directory.
    - image_dimension: Dimensions of the image (default is [2048, 2048]).
    """

    # Import necessary libraries
    from os import listdir
    from os.path import isfile, join
    import tifffile
    from xml.dom import minidom
    import pandas as pd
    import numpy as np
    import os
    from tifffile import imread
    from tqdm import tqdm
    import re
    import shutil

    # Refactor input directories for compatibility (especially with Linux)
    refactored_dirs = [dir_path.replace("%20", " ") for dir_path in input_dirs]

    # Iterate through each input directory
    for idx, dir_path in enumerate(refactored_dirs):

        # Get list of files in the directory
        files = os.listdir(dir_path)
        
        # Filter for TIFF files that are not deconvolved
        tif_files = [file for file in files if 'dw' not in file and '.tif' in file and '.txt' not in file]
        
        # Split filenames to get the regions
        split_underscore = pd.DataFrame(tif_files)[0].str.split('--', expand=True)
        unique_regions = list(split_underscore[0].unique())

        # If the scan is large, it may be divided into multiple regions
        for region in unique_regions:
            region_tif_files = [file for file in tif_files if region in file]
            base_index = str(idx + 1)
            split_underscore = pd.DataFrame(region_tif_files)[0].str.split('--', expand=True)
            
            # Extract tiles information
            tiles = sorted(split_underscore[1].unique())
            tiles_df = pd.DataFrame(tiles)
            tiles_df['indexNumber'] = [int(tile.split('e')[-1]) for tile in tiles_df[0]]
            tiles_df.sort_values(by=['indexNumber'], ascending=True, inplace=True)
            tiles_df.drop('indexNumber', 1, inplace=True)
            tiles = list(tiles_df[0])
            
            # Extract channels information
            channels = split_underscore[3].unique()
            
            # Determine the output directory based on the region
            if len(unique_regions) == 1:
                output_dir = output_dir_prefix
            else:
                output_dir = f"{output_dir_prefix}_R{region.split('Region')[1].split('_')[0]}"
            mipped_output_dir = f"{output_dir}/preprocessing/mipped/"
            
            # Create directory if it doesn't exist
            if not os.path.exists(mipped_output_dir):
                os.makedirs(mipped_output_dir)

            for base_idx, base in enumerate(sorted(base_index)):
                if not os.path.exists(f"{mipped_output_dir}/Base_{base}"):
                    os.makedirs(f"{mipped_output_dir}/Base_{base}")
                try:
                    metadata_file = join(dir_path, 'Metadata', [file for file in os.listdir(join(dir_path, 'Metadata')) if region in file][0])
                    if not os.path.exists(join(mipped_output_dir, f"Base_{base}", 'MetaData')):
                        os.makedirs(join(mipped_output_dir, f"Base_{base}", 'MetaData'))
                    shutil.copy(metadata_file, join(mipped_output_dir, f"Base_{base}", 'MetaData'))
                except FileExistsError:
                    pass

                # Maximum Intensity Projection (MIP) for each tile
                for _tile in tqdm(range(len(tiles))):
                    tile = tiles[_tile]
                    tile_for_name = re.split('(\d+)', tile)[1]
                    existing_files = [file for file in os.listdir(f"{mipped_output_dir}/Base_{base}") if str(tile_for_name) in file]
                    
                    # Ensure that we don't overwrite existing files
                    if len(existing_files) < len(channels):
                        tile_tif_files = [file for file in region_tif_files if f"{tile}--" in file]
                        for channel_idx, channel in enumerate(sorted(list(channels))):
                            channel_tif_files = [file for file in tile_tif_files if str(channel) in file]
                            max_intensity = np.zeros(image_dimension)
                            for file in channel_tif_files:
                                try:
                                    im_array = imread(f"{dir_path}/{file}")
                                except:
                                    print('Image corrupted, reading black file instead.')
                                    im_array = np.zeros(image_dimension)
                                max_intensity = np.maximum(max_intensity, im_array)
                            max_intensity = max_intensity.astype('uint16')
                            tifffile.imwrite(f"{mipped_output_dir}/Base_{base}/Base_{base}_s{tile_for_name}_{channel}", max_intensity)







def leica_OME_tiff(directory_base, output_directory):
    """
    Convert Leica TIFF files to OME-TIFF format.
    
    Args:
    - directory_base: Base directory containing the TIFF files.
    - output_directory: Directory to save the converted OME-TIFF files.
    
    Returns:
    None. Writes the OME-TIFF images to the designated output directory.
    """
    
    import tifffile
    import numpy as np
    import os
    from os.path import join
    import tifffile
    import os
    from os import listdir
    import pandas as pd
    import numpy as np
    from xml.dom import minidom
    from pathlib import Path
    from tqdm import tqdm

    folders = os.listdir(directory_base)
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    for folder in folders:
        exported_directory = join(directory_base,folder)
        onlyfiles = listdir(exported_directory)
        onlytifs =  [k for k in onlyfiles if '.tif' in k]
        onlyfiles_df = pd.DataFrame(onlytifs)
        onlyfiles_split_tiles = onlyfiles_df[0].str.split('_s',expand=True)
        onlyfiles_split_channel = onlyfiles_split_tiles[1].str.split('_',expand=True)

        tiles = list(np.unique(onlyfiles_split_tiles[1].str.split('_',expand=True)[0]))
        tiles_df=pd.DataFrame(tiles)
        tiles_df['indexNumber'] = [int(i.split('e')[-1]) for i in tiles_df[0]]
        # Perform sort of the rows
        tiles_df.sort_values(by = ['indexNumber'], ascending = [True], inplace = True)
        # Deletion of the added column
        tiles_df.drop('indexNumber', 1, inplace = True)
        tiles = list(tiles_df[0])
        channels = list(np.unique(onlyfiles_split_channel[1]))
        rounds = list(np.unique(onlyfiles_split_tiles[0]))
        
        
        metadatafiles = listdir(join(exported_directory, 'MetaData'))
        metadatafiles =  [k for k in metadatafiles if 'IOManagerConfiguation.xlif' not in k]

        for p, meta in enumerate(metadatafiles):
            print(meta)
            mydoc = minidom.parse(join(exported_directory, 'MetaData',meta) )
            tile =[]
            x =[]
            y =[]
            items = mydoc.getElementsByTagName('Tile')
            for el, elem in enumerate(items):
                tile.append(el)
                x.append(float(elem.attributes['PosX'].value))
                y.append(float(elem.attributes['PosY'].value))
            unique_tiles = list(np.unique(tile))
            x_reformatted = (x[:len(unique_tiles)])    
            y_reformatted = (y[:len(unique_tiles)])     
            dictionary = {'x': x_reformatted, 'y': y_reformatted}  

            df = pd.DataFrame(dictionary)
            df['x'] =((df.x-np.min(df.x))/.000000321) + 1
            df['y'] =((df.y-np.min(df.y))/.000000321) + 1
            positions = np.array(df).astype(int)
            df.to_csv(directory_base +'/'+ folder + '/coords.csv')
            
        with tifffile.TiffWriter((output_directory +'/'+ folder + '.ome.tiff'), bigtiff=True) as tif:
            for i in tqdm(range(len(tiles))):
                position = positions[i]
                tile = tiles[i]

                tile_filtered = [k for k in onlytifs if 's'+tile+'_' in k]
                tile_filtered =  [k for k in tile_filtered if '._' not in k]

                stacked = np.empty((5, 2048, 2048))
                for n,image_file in enumerate(sorted(tile_filtered)):
                    try: 
                        image_int = tifffile.imread(join(exported_directory,image_file))
                    except IndexError: 
                        image_int = np.empty((2048, 2048))
                    stacked[n] = image_int.astype('uint16')
                pixel_size = 0.1625
                metadata = {
                                'Pixels': {
                                    'PhysicalSizeX': pixel_size,
                                    'PhysicalSizeXUnit': 'µm',
                                    'PhysicalSizeY': pixel_size,
                                    'PhysicalSizeYUnit': 'µm'
                                },
                                'Plane': {
                                    'PositionX': [position[0]*pixel_size]*stacked.shape[0],
                                    'PositionY': [position[1]*pixel_size]*stacked.shape[0]
                                }

                            }
                tif.write(stacked.astype('uint16'),metadata=metadata)


def ashlar_wrapper(
    files, 
    output='', 
    align_channel=1, 
    flip_x=False, 
    flip_y=True, 
    output_channels=None, 
    maximum_shift=500, 
    filter_sigma=5.0, 
    filename_format='Round{cycle}_{channel}.tif',
    pyramid=False,
    tile_size=None,
    ffp=False,
    dfp=False,
    plates=False,
    quiet=False,
    version=False):
    """
    Wrapper function for the Ashlar tool for image alignment and mosaicking.
    
    Args:
    - files (list): List of input image file paths.
    - output (str): Directory path to save the output images.
    - ... (other args): Various configuration parameters for image processing.
    
    Returns:
    - int: 1 for errors, else the output of the Ashlar processing function.
    """
    
    ashlar.configure_terminal()
    
    filepaths = files
    output_path = pathlib.Path(output)

    import warnings
    warnings.filterwarnings("ignore")   

    # make directory
    if not os.path.exists(output):
        os.makedirs(output)

    if tile_size and not pyramid:
        ashlar.print_error("--tile-size can only be used with --pyramid")
        return 1
    if tile_size is None:
        # Implement default value logic as mentioned in argparser setup above.
        tile_size = tile_size

    ffp_paths = ffp
    if ffp_paths:
        if len(ffp_paths) not in (0, 1, len(filepaths)):
            ashlar.print_error(
                "Wrong number of flat-field profiles. Must be 1, or {}"
                " (number of input files)".format(len(filepaths))
            )
            return 1
        if len(ffp_paths) == 1:
            ffp_paths = ffp_paths * len(filepaths)

    dfp_paths = dfp
    if dfp_paths:
        if len(dfp_paths) not in (0, 1, len(filepaths)):
            ashlar.print_error(
                "Wrong number of dark-field profiles. Must be 1, or {}"
                " (number of input files)".format(len(filepaths))
            )
            return 1
        if len(dfp_paths) == 1:
            dfp_paths = dfp_paths * len(filepaths)

    aligner_args = {}
    aligner_args['channel'] = align_channel
    aligner_args['verbose'] = not quiet
    aligner_args['max_shift'] = maximum_shift
    aligner_args['filter_sigma'] = filter_sigma

    mosaic_args = {}
    if output_channels:
        mosaic_args['channels'] = output_channels
    if pyramid:
        mosaic_args['tile_size'] = tile_size
    if quiet is False:
        mosaic_args['verbose'] = True

    try:
        if plates:
            return ashlar.process_plates(
                filepaths, output_path, filename_format, flip_x,
                flip_y, ffp_paths, dfp_paths, aligner_args, mosaic_args,
                pyramid, quiet
            )
        else:
            mosaic_path_format = str(output_path / filename_format)
            return ashlar.process_single(
                filepaths, mosaic_path_format, flip_x, flip_y,
                ffp_paths, dfp_paths, aligner_args, mosaic_args, pyramid,
                quiet
            )
    except ashlar.ProcessingError as e:
        ashlar.print_error(str(e))
        return 1



def reshape_split(image: np.ndarray, kernel_size: tuple):
    """
    Reshape the input image into smaller tiles of specified size.
    
    Args:
    - image (np.ndarray): 2D array representing the image.
    - tile_size (tuple): Desired dimensions for the tiles (height, width).
    
    Returns:
    - np.ndarray: 4D array representing tiles of the image.
    """    
    img_height, img_width = image.shape
    tile_height, tile_width = kernel_size
    
    tiled_array = image.reshape(img_height // tile_height, 
                               tile_height, 
                               img_width // tile_width, 
                               tile_width)
    
    tiled_array = tiled_array.swapaxes(1,2)
    return tiled_array

    


def tile_stitched_images(image_path,outpath, tile_dim=2000, file_type = 'tif', old_stiched_name = False):

    """
    Tiles stitched images from a directory and saves them with a specific naming convention.
    
    Args:
    - image_directory (str): Directory containing the stitched images.
    - output_directory (str): Directory to save the tiled images.
    - tile_dim (int): Dimension for tiling. Default is 2000.
    - file_type (str): Type of the image file. Default is 'tif'.
    - old_stitched_naming (bool): Flag to handle old naming convention. Default is False.
    """
    
    if not os.path.exists(outpath):
            os.makedirs(outpath)
            
    images = os.listdir(image_path)
    images =  [k for k in images if '._' not in k]
    
    if file_type=='mat':
        images =  [k for k in images if '.tif.mat' in k] 
    else: 
        images =  [k for k in images if '.tif' in k] 

    for image_file in sorted(images):
        try: 
            if file_type == 'mat':
                image = mat73.loadmat(image_path +'/'+ image_file)['I']
                cycle = ''.join(filter(str.isdigit, image_file.split('_')[1]))
                channel = ''.join(filter(str.isdigit, image_file.split('_')[2].split('-')[1].split('.')[0]))
            else:
                if old_stiched_name == True:
                    print('old names')
                    image = tifffile.imread(image_path +'/'+ image_file)
                    cycle = str(int(''.join(filter(str.isdigit, image_file.split('_')[1])))-1)
                    channel = str(int(''.join(filter(str.isdigit, image_file.split('-')[1])))-1)
                    print(cycle)
                    print(channel)
                else: 
                    image = tifffile.imread(image_path +'/'+ image_file)
                    cycle = ''.join(filter(str.isdigit, image_file.split('_')[0]))
                    channel = ''.join(filter(str.isdigit, image_file.split('_')[1]))

           
            
            print('tiling: ' + image_file)
            
            image_pad = cv2.copyMakeBorder( image, top = 0, bottom =math.ceil(image.shape[0]/tile_dim)*tile_dim-image.shape[0], left =0, right = math.ceil(image.shape[1]/tile_dim)*tile_dim-image.shape[1], borderType = cv2.BORDER_CONSTANT)
            image_split = reshape_split(image_pad,(tile_dim,tile_dim))
            nrows, ncols, dim1, dim2 = image_split.shape
            x = []
            y = []
            directory = outpath +'/'+'Base_'+str(int(cycle)+1)+'_stitched-'+str(int(channel)+1) 
            if not os.path.exists(directory):
                os.makedirs(directory) 
            count = 0
            for i in range(nrows):
                for j in range(ncols):
                    count = count+1                
                    x.append(j*tile_dim)
                    y.append(i*tile_dim)
                    
                    tifffile.imwrite(directory + '/' +'tile'+str(count)+'.tif',image_split[i][j])
        except KeyError:
            continue
                
    tile_pos = pd.DataFrame()
    tile_pos['x'] = x
    tile_pos['y'] = y

    tile_pos.to_csv(outpath+'/'+'tilepos.csv', header=False, index=False)
    return
    
    
def preprocessing_main_leica(input_dirs, 
                            output_location,
                            regions_to_process = 2, 
                            align_channel = 4, 
                            tile_dimension = 6000, 
                            mip = True):
    """
    Main function to preprocess Leica microscopy images.

    Args:
    - input_dirs (str): Directories containing input images.
    - output_location (str): Base output directory.
    - regions_to_process (int): Number of regions to process. Default is 2.
    - align_channel (int): Channel to use for alignment. Default is 4.
    - tile_dimension (int): Dimension for tiling. Default is 6000.
    - mip (bool): Flag to perform maximum intensity projection. Default is True. Use false for pre-mipped images
    """
    
    # Maximum Intensity Projection
    if mip == True:
        leica_mipping(input_dirs=input_dirs, output_dir_prefix = output_location)
    else: 
        print('not mipping')
        
    if regions_to_process > 1:
        for i in range(regions_to_process):
            path = output_location +'_R'+str(i+1)
            
            # create leica OME_tiffs
            leica_OME_tiff(directory_base = path+'/preprocessing/mipped/', 
                                            output_directory = path+'/preprocessing/OME_tiffs/')
            
            # align and stitch images
            OME_tiffs = os.listdir(path+'/preprocessing/OME_tiffs/')
            OME_tiffs = [path+'/preprocessing/OME_tiffs/' + sub for sub in OME_tiffs]
            ashlar_wrapper(files = OME_tiffs, 
                                            output = path+'/preprocessing/stitched/', 
                                            align_channel=align_channel)
            
            # retile stitched images
            tile_stitched_images(image_path = path+'/preprocessing/stitched/',
                                    outpath = path+'/preprocessing/ReslicedTiles/', 
                                    tile_dim=tile_dimension)

    
    else: 
        path = output_location

        # create leica OME_tiffs
        leica_OME_tiff(directory_base = path+'/preprocessing/mipped/', 
                                        output_directory = path+'/preprocessing/OME_tiffs/')

        # align and stitch images
        OME_tiffs = os.listdir(path+'/preprocessing/OME_tiffs/')
        OME_tiffs = [path+'/preprocessing/OME_tiffs/' + sub for sub in OME_tiffs]

        ashlar_wrapper(files = OME_tiffs, 
                                        output = path+'/preprocessing/stitched/', 
                                        align_channel=align_channel)

        # retile stitched images
        tile_stitched_images(image_path = path+'/preprocessing/stitched/',
                                outpath = path+'/preprocessing/ReslicedTiles/', 
                                tile_dim=tile_dimension)
    return






def process_czi(input_file, outpath, mip=True, cycle=0, tile_size_x=2048, tile_size_y=2048):
    """
    Process CZI files, apply maximum intensity projection (if specified), 
    and create an associated XML with metadata.
    
    Parameters:
    - input_file: Path to the input CZI file.
    - outpath: Directory where the processed images and XML will be saved.
    - mip: Boolean to decide whether to apply maximum intensity projection. Default is True.
    - cycle: Int to specify the cycle number. Default is 0.
    - tile_size_x: Size of the tile in X dimension. Default is 2048.
    - tile_size_y: Size of the tile in Y dimension. Default is 2048.
    
    Returns:
    - A string indicating that processing is complete.
    """
    
    # import packages 
    import os
    import xml.etree.ElementTree as ET
    from aicspylibczi import CziFile
    import aicspylibczi
    from xml.dom import minidom
    import numpy as np
    from tqdm import tqdm
    import pandas as pd
    import tifffile
    
    # Create the output directory if it doesn't exist.
    if not os.path.exists(outpath):
        os.makedirs(outpath)


    # Load the CZI file and retrieve its dimensions.
    czi = aicspylibczi.CziFile(input_file)
    dimensions = czi.get_dims_shape() 
    chsize = dimensions[0]['C'][1]
    try:
        msize=dimensions[0]['M'][1]
    except:
        msize=0 
    ssize = dimensions[0]['S'][1]

    # Check if mip is True and cycle is not zero.
    if mip and cycle != 0:
        # Initialize placeholders for metadata.
        Bxcoord = []
        Bycoord = []
        Btile_index = []
        filenamesxml = []
        Bchindex = []

        # Loop through each mosaic tile and each channel.
        for m in tqdm(range(0, msize)):
            for ch in range (0, chsize):
                # Get metadata and image data for the current tile and channel.
                meta = czi.get_mosaic_tile_bounding_box(M=m, Z=0, C=ch)
                img, shp = czi.read_image(M=m, C=ch)
                
                # Apply maximum intensity projection.
                IM_MAX = np.max(img, axis=3)
                IM_MAX = np.squeeze(IM_MAX, axis=(0,1,2,3))
                
                # Construct filename for the processed image.
                n = str(0)+str(m+1) if m < 9 else str(m+1)
                filename = 'Base_' + str(cycle) + '_c' + str(ch+1) + 'm' + str(n) + '_ORG.tif'
                
                # Save the processed image.
                
                tifffile.imwrite(outpath + filename, IM_MAX.astype('uint16'))
                
                # Append metadata to the placeholders.
                Bchindex.append(ch)
                Bxcoord.append(meta.x)
                Bycoord.append(meta.y)
                Btile_index.append(m)
                filenamesxml.append(filename)

        # Adjust the XY coordinates to be relative.
        nBxcord = [x - min(Bxcoord) for x in Bxcoord]
        nBycord = [y - min(Bycoord) for y in Bycoord]
        
        # Create a DataFrame to organize the collected metadata.
        metadatalist = pd.DataFrame({
            'Btile_index': Btile_index, 
            'Bxcoord': nBxcord, 
            'Bycoord': nBycord, 
            'filenamesxml': filenamesxml,
            'channelindex': Bchindex
        })
        
        metadatalist = metadatalist.sort_values(by=['channelindex','Btile_index'])
        metadatalist.reset_index(drop=True)

        # Initialize the XML document structure.
        export_doc = ET.Element('ExportDocument')
        
        # Populate the XML document with metadata.
        for index, row in metadatalist.iterrows():
            image_elem = ET.SubElement(export_doc, 'Image')
            filename_elem = ET.SubElement(image_elem, 'Filename')
            filename_elem.text = row['filenamesxml']
            
            bounds_elem = ET.SubElement(image_elem, 'Bounds')
            bounds_elem.set('StartX', str(row['Bxcoord']))
            bounds_elem.set('SizeX', str(tile_size_x))
            bounds_elem.set('StartY', str(row['Bycoord']))
            bounds_elem.set('SizeY', str(tile_size_y))
            bounds_elem.set('StartZ', '0')
            bounds_elem.set('StartC', '0')
            bounds_elem.set('StartM', str(row['Btile_index']))
            
            zoom_elem = ET.SubElement(image_elem, 'Zoom')
            zoom_elem.text = '1'

        
        # Save the constructed XML document to a file.
        xml_str = ET.tostring(export_doc)
        with open(outpath + 'Base_' + str(cycle) + '_info.xml', 'wb') as f:
            f.write(xml_str)

    return "Processing complete."
