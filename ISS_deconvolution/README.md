In the `ISS_deconvolution` module, we apply denoising/deconvolution methods to sharpen the images of the rolling circle products, before proceeding to decoding. This is especially advantageous in case of low magnification imaging and/or in signal-crowded tissues, because deconvolution helps separating nearby dots and their respective signals. 

Deconvolution is implemented with the `flowdec` library ([https://github.com/hammerlab/flowdec]()) for tensorflow-accelerated image deconvolution. The library has decent quality documentation, and we suggest you refer to it for advanced tweaking and options. Nevertheless, the library is not currently maintained.

Our plan is, in the long run, to move away from `flowdec`  and use something else, but we haven't made up our mind about this yet.




Installation: You can download installation.sh and just execute it. This will install the following modules:

`ISS_preprocessing`

`ISS_decoding`

`ISS_postprocessing`

`ISS_deconvolution`

`ISS_CARE`

`ISS_RNA_probedesign`

Please refer to the respective `readme.md` files for further instructions.