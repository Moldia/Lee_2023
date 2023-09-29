from __future__ import print_function, unicode_literals, absolute_import, division
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
%config InlineBackend.figure_format = 'retina'
import tifffile
import skimage

import os
from tifffile import imread, imsave
from csbdeep.utils import Path, download_and_extract_zip_file, plot_some
from csbdeep.io import save_tiff_imagej_compatible
from csbdeep.models import CARE
import shutil

def ISS_CARE(directory,
             basedir_model,
             model_name,
             output_dir,
             DAPI_ch=5
            ):
    
    """
    directory  = type: str. this is the path to your /ReslicedTiles/ 
    folder containing the raw mipped images to deconvolve.

    basedir_model = type: str. The path to the folder containing the pre-trained models. 
    
    model_name = type: str. The name of the model to apply. 

    output_dir = type: str. This is the output directory for the denoised images. 
    A /preprocessing/ReslicedTiles/ folder will be created here, along with specific subfolders 
    reflecting the folder tree of the input path.

    DAPI_ch = type: int. Number of the channel containing the DAPI images. 
    This is to exclude DAPI from the CARE processing.
    """
    
    axes = 'YX'
    directoriestopredict = []
    dapifolders=[]
    for root, subdirectories, files in os.walk(directory):
        for subdirectory in subdirectories:
            if subdirectory[-1]==str(DAPI_ch):
                dapifolders.append(subdirectory)
                #print ('DAPI')
            else:
                #print ('spots')
                directoriestopredict.append(subdirectory)

            #filestopredict = []
            #for file in files:
            #    if file.endswith(".tif"):
            #        filestopredict.append(file)
        #    print(os.path.join(root, file))
    filestopredict = []
    for file in os.listdir(directory+subdirectory):
        if file.endswith(".tif"):
            filestopredict.append(file)

    dapitopredict = []
    for file in os.listdir(directory+subdirectory):
        if file.endswith(".tif"):
            dapitopredict.append(file)
    
    model = CARE(config=None, name=model_name, basedir=basedir_model)
    
    output_dir = output_dir + model.name + '/ReslicedTiles/'
    if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
    for fold in directoriestopredict:
        print (fold)
        for file in filestopredict:
            print (file)
            #print (fold+file)
            x = imread(directory+'/'+fold+'/'+file)
        #restored = model.predict(x, axes))
            restored = model.predict(x, axes, n_tiles=(1,2))
            if not os.path.exists(output_dir+fold):
                os.makedirs(output_dir+fold)
            #restored=skimage.img_as_float32(restored)
            #skimage.io.imsave(output_dir + fold +'/'+ file, restored)
            tifffile.imsave(output_dir + fold +'/'+ file, restored.astype('uint16'))
    #restored = model.predict(x, axes)    
    
    #Copy the DAPI folders without CAREing them
    for dapi in dapifolders:
        try:
            dapisource=directory+dapi
            dapidest=output_dir+dapi
            shutil.copytree(dapisource, dapidest)
        except:
            continue

    # Copy the metadata into the CAREd folders
    src_metadata = directory+'tilepos.csv'
    dst_metadata = output_dir+'tilepos.csv'
    shutil.copyfile(src_metadata, dst_metadata)    