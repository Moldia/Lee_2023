from os import listdir
from os.path import isfile, join
from xml.dom import minidom
import pandas as pd
import numpy as np
import re
import shutil
import os
import tifffile
import matplotlib.pyplot as plt
from flowdec.nb import utils as nbutils 
from flowdec import psf as fd_psf
from flowdec import data as fd_data
from scipy import ndimage
import dask
import dask.array as da
import tensorflow as tf
from flowdec.restoration import RichardsonLucyDeconvolver
from skimage import io
from pathlib import Path
import operator
from flowdec import data as fd_data
from tqdm import tqdm
import os
import xml.etree.ElementTree as ET
from aicspylibczi import CziFile


'''
#THIS IS AN EXAMPLE OF PSF DATA
PSF_metadata = {'na':0.8,
'm':20,
'ni0':1.42,
'res_lateral':0.419,
'res_axial':1.718,
 'channels':{
 'AF750':{
    'wavelength':.773},
  'Cy5':{
    'wavelength':.673},
  'Cy3':{
    'wavelength':.561},
  'AF488':{
    'wavelength':.519},
 'DAPI':{
    'wavelength':.465}
 }
}
m= magnification
ni0 = refraction coefficient of the immersion medium
res_lateral = resolution in xy
res_axial = resolution in z
'''








def deconvolve_czi(input_file, outpath, image_dimensions=[2048, 2048], PSF_metadata=None, chunk_size=None,  mip=True, cycle=0, tile_size_x=2048, tile_size_y=2048):

    """
    Process CZI files, deconvolve the image stacks, apply maximum intensity projection (if specified), 
    and create an associated XML with metadata.
    
    Parameters:
    - input_file: Path to the input CZI file.
    - outpath: Directory where the processed images and XML will be saved.
    - chunk_size= [x,y] where x and y are the size of the chunks the image needs to be cut into for small GPU processing
    - mip: Boolean to decide whether to apply maximum intensity projection. Default is True.
    - cycle: Int to specify the cycle number. Default is 0.
    - tile_size_x: Size of the tile in X dimension. Default is 2048.
    - tile_size_y: Size of the tile in Y dimension. Default is 2048.
    
    Returns:
    - A string indicating that processing is complete.
    """
	
def customcopy(src, dst):
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    shutil.copyfile(src, dst)

    def cropND(img, bounding):
        start = tuple(map(lambda a, da: a//2-da//2, img.shape, bounding))
        end = tuple(map(lambda a, da: a+da, start, bounding))
        slices = tuple(map(slice, start, end))
        return img[slices]

    # Load the CZI file and retrieve its dimensions.
    czi = aicspylibczi.CziFile(input_file)
    dimensions = czi.get_dims_shape() 
    chsize = dimensions[0]['C'][1]
    msize = dimensions[0]['M'][1]
    z_size=dimensions[0]['Z'][1]

# Create the output directory if it doesn't exist.
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    # Check if mip is and cycle is not zero.
    if cycle != 0:
        Bxcoord = []
        Bycoord = []
        Btile_index = []
        filenamesxml = []
        Bchindex = []
        
        psf_dict = {}
        for idx, channel in enumerate(sorted(PSF_metadata['channels'])):
            psf_dict[idx] = fd_psf.GibsonLanni(
                na=float(PSF_metadata['na']),
                m=float(PSF_metadata['m']),
                ni0=float(PSF_metadata['ni0']),
                res_lateral=float(PSF_metadata['res_lateral']),
                res_axial=float(PSF_metadata['res_axial']),
                wavelength=float(PSF_metadata['channels'][channel]['wavelength']),                
                size_x=tile_size_x,
                size_y=tile_size_y,
                size_z=z_size  # Use the Z dimension from the CZI file
            ).generate()

        # Process each channel and mosaic tile
        for m in tqdm(range(0, msize)):
            for ch in range(0, chsize):
                # Read image data for the current mosaic tile and channel
                img, shp = czi.read_image(M=m, C=ch)
                img=np.squeeze(img, axis=(0,1,2,4))
                
                # Read metadata for the current mosaic tile and channel

                meta = czi.get_mosaic_tile_bounding_box(M=m, Z=0, C=ch)
           

                # Instantiate the Richardson-Lucy deconvolution algorithm
                algo = RichardsonLucyDeconvolver(img.ndim, pad_mode="2357", pad_min=(6,6,6))

                # Check if chunk_size is provided
                if chunk_size:
                    # Define chunk dimensions
                    chunked_dims = (z_size, chunk_size[0], chunk_size[1])

                    # Convert image data to a chunked dask array
                    arr = da.from_array(img, chunks=chunked_dims)
                    cropped_kernel = cropND(psf_dict[ch], chunked_dims)

                    # Define deconvolution function for chunks
                    def deconv(chunk):
                        tmp = algo.initialize().run(fd_data.Acquisition(data=chunk, kernel=cropped_kernel), 50)
                        return tmp.data 

                    # Apply chunked deconvolution
                    deconvolved = arr.map_overlap(deconv, depth=(6,6,6), boundary='reflect', dtype='uint16').compute(num_workers=1)
                else:
                    # Regular deconvolution for the entire image
                    deconvolved = algo.initialize().run(fd_data.Acquisition(data=img, kernel=psf_dict[ch]), 50)
                if chsize != len(PSF_metadata['channels']):
                    raise ValueError("Mismatch between CZI file channels and PSF_metadata channels.")
                #print(deconvolved.data.shape)
                # Check if mip (max intensity projection) is enabled
                if mip:
                    processed_img = np.max(deconvolved.data, axis=0).astype('uint16')
                else:
                    processed_img = deconvolved.data.astype('uint16')
            
                # Construct filename for the processed image
                n = str(0)+str(m+1) if m < 9 else str(m+1)
                filename = f'Base_{cycle}_c{ch+1}m{n}_ORG.tif'

                # Save the processed image
                tifffile.imwrite(os.path.join(outpath, filename), processed_img)
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
    




def deconvolve_leica(input_dirs, output_dir_prefix, image_dimensions=[2048, 2048], PSF_metadata=None, chunk_size=None, mip=True):
    """
    Process the images from the given directories.

    Parameters:
    - input_dirs: List of directories containing the images to process.
    - output_dir_prefix: Prefix for the output directories.
    - image_dimensions: Dimensions of the images (default: [2048, 2048]).
    - PSF_metadata: Metadata for Point Spread Function (PSF) generation.
    - chunk_size [x,y]: Size of chunks for processing. If None, the entire image is processed.
                  Small GPUs will require chunking. Enable if you run out of gRAM
    - mip: Boolean to decide whether to apply maximum intensity projection. Default is True. 
           If mip=false the stack is deconvolved but saved as an image stack without projecting it

    Returns:
    None. Processed images are saved in the output directories.
    """
    def cropND(img, bounding):
        start = tuple(map(lambda a, da: a//2-da//2, img.shape, bounding))
        end = tuple(map(lambda a, da: a+da, start, bounding))
        slices = tuple(map(slice, start, end))
        return img[slices]
    
    if PSF_metadata is None:
        raise ValueError("PSF_metadata is required to generate PSF.")

    def custom_copy(src, dest):
        """Custom function to copy a file to a destination."""
        if os.path.isdir(dest):
            dest = os.path.join(dest, os.path.basename(src))
        shutil.copyfile(src, dest)

    # Preprocess directory names (replace placeholders)
    processed_dirs = [dir_name.replace("%20", " ") for dir_name in input_dirs]  # Replace placeholder for spaces

    # Loop through directories to process images
    for dir_index, dir_path in enumerate(processed_dirs):
        tif_files = [f for f in os.listdir(dir_path) if '.tif' in f and 'dw' not in f and '.txt' not in f]
        
        # Extract unique region identifiers from filenames
        regions = pd.DataFrame(tif_files)[0].str.split('--', expand=True)[0].unique()

       
        # Process each region
        for region in regions:
            print('Counting the regions and organizing the files')
            filtered_tifs = [f for f in tif_files if region in f]
            base_num = str(dir_index + 1)
            
            # Extract unique tile identifiers
            tiles_df = pd.DataFrame(filtered_tifs)[0].str.split('--', expand=True)[1]
            tiles_df = tiles_df.str.extract('(\d+)')[0].sort_values().unique()

            # Determine output directory based on number of regions
            if len(regions) == 1:
                output_directory = output_dir_prefix
                mipped_directory = os.path.join(output_directory, 'preprocessing', 'mipped')
            else:
                output_directory = f"{output_dir_prefix}_R{region.split('Region')[1].split('_')[0]}"
                mipped_directory = os.path.join(output_directory, 'preprocessing', 'mipped')

            # Ensure the output directory exists
            os.makedirs(mipped_directory, exist_ok=True)
            
            # Process each cycle
            for base in sorted(base_num):
                base_directory = os.path.join(mipped_directory, f'Base_{base}')
                os.makedirs(base_directory, exist_ok=True)
                
                # Copy metadata if it exists
                print ('Extracting metadata')
                metadata_dir = os.path.join(dir_path, 'Metadata')
                metadata_files = [f for f in os.listdir(metadata_dir) if region in f]
                if metadata_files:
                    os.makedirs(os.path.join(base_directory, 'MetaData'), exist_ok=True)
                    custom_copy(os.path.join(metadata_dir, metadata_files[0]), os.path.join(base_directory, 'MetaData'))
                
                # Extracts the first tile to calculate its Z depth
                print ('Calculating the PSF')
                sample_tile=[f for f in filtered_tifs if f"--Stage00--" in f]
                size_z = int(len(sample_tile) / len(PSF_metadata['channels']))
                # Generate PSFs for each channel outside the tile loop
                psf_dict = {}
                for channel in PSF_metadata['channels']:
                    psf_dict[channel] = fd_psf.GibsonLanni(
                        na=float(PSF_metadata['na']),
                        m=float(PSF_metadata['m']),
                        ni0=float(PSF_metadata['ni0']),
                        res_lateral=float(PSF_metadata['res_lateral']),
                        res_axial=float(PSF_metadata['res_axial']),
                        wavelength=float(PSF_metadata['channels'][channel]['wavelength']),
                        size_x=image_dimensions[0],
                        size_y=image_dimensions[1],
                        size_z=size_z
                    ).generate()
                # Process each tile within the base
                print ('Deconvolving Cycle: '+ base)
                for tile in tqdm(sorted(tiles_df, key=int)):
                #for tile in sorted(tiles_df, key=int):
                    tile_files = [f for f in filtered_tifs if f"--Stage{tile}--" in f]
                    for channel in sorted(PSF_metadata['channels']):
                        output_file_path = os.path.join(base_directory, f'Base_{base}_s{tile}_C0{channel}.tif')

                        if os.path.exists(output_file_path):
                            print(f"File {output_file_path} already exists. Skipping this tile for channel {channel}.")
                            continue
                        
                        channel_files = [f for f in tile_files if f"--C0{channel}" in f]
                        stacked_images = np.stack([tifffile.imread(os.path.join(dir_path, f)) for f in channel_files])

                        if chunk_size:
                            # Convert image data to a chunked dask array
                            # Set chunk size to include z dimension, according to the specified XY
                            adjusted_chunk_size = (size_z, chunk_size[0], chunk_size[1])

                            arr = da.from_array(stacked_images, chunks=adjusted_chunk_size)
                            cropped_kernel = cropND(psf_dict[channel], adjusted_chunk_size)
                            algo = RichardsonLucyDeconvolver(stacked_images.ndim, pad_mode="2357", pad_min=(6,6,6))

                            def deconv(chunk):
                                tmp = algo.initialize().run(fd_data.Acquisition(data=chunk, kernel=cropped_kernel), 50)
                                return tmp.data 

                            deconvolved = arr.map_overlap(deconv, depth=(6,6,6), boundary='reflect', dtype='uint16').compute(num_workers=1)
                            #result_maxproj = np.max(result_overlap, axis=0).astype('uint16')
                        else:
                            # Regular deconvolution for the entire image
                            algo = RichardsonLucyDeconvolver(stacked_images.ndim, pad_mode="2357", pad_min=(6,6,6))
                            deconvolved = algo.initialize().run(fd_data.Acquisition(data=stacked_images, kernel=psf_dict[channel]), 50)
                            #result_maxproj = np.max(deconvolved.data, axis=0).astype('uint16')
                        
                        if mip:
                            processed_img = np.max(deconvolved.data, axis=0).astype('uint16')
                        else:
                            processed_img = np.asarray(deconvolved.data).astype('uint16')

                        tifffile.imwrite(os.path.join(base_directory, f'Base_{base}_s{tile}_C0{channel}.tif'), processed_img)

    return None
