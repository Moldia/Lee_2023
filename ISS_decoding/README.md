The `ISS_decoding`  we make use of the starfish library to extract the information contained in our images. starfish is a Python library for processing images of image-based spatial transcriptomics. You can read more about it here: https://spacetx-starfish.readthedocs.io/en/latest/

The module contains function to first format our preprocessed images to a starfish-compatible format (SpaceTx).
In the next steps we then proceed to the actual decoding of the data from the SpaceTx images, and a further set of functions allows the user to plot the decoding results and filter the data appropriately.


