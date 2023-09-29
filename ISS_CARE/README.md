# ISS_CARE
Content Aware Image Restoration (CARE) is a ML-based method for image denoising. When appropriately trained, it can be used as a very fast alternative to image deconvolution. Its main advantage is that it can work directly on projected images, significantly reducing the computing requirements. This speeds up the denoising process by a factor of at least a 100 compared to Flowdec (which is already a faster method than most available), although we acknowledge it is likely less accurate than "true" deconvolution. From our benchmarking, CARE seems to be the way to go for most experiments. We advise the user to do a test run on a small sample with both Flowdec and CARE, compare the results and proceed using CARE only if the results are satisfyingly  similar.

You can find the examples on how to use it in the folder "Notebooks".

You can find the pretrained models in the folder "models"
