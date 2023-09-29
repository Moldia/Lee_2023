from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.23'
DESCRIPTION = 'Package used to deconvolve raw ISS images from Leica and Zeiss microscopes using GPU'
LONG_DESCRIPTION = 'This package uses flowdec to deconvolve and project images for downstream processing. Replaces the mipping function from preprocessing'

# Setting up
setup(
    name="ISS_deconvolution",
    version=VERSION,
    author="Marco Grillo",
    author_email="<marco.grillo@scilifelab.se>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=[],
    keywords=['python', 'spatial transcriptomics', 
            'spatial resolved transcriptomics', 
            'in situ sequencing', 
            'ISS','decoding'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Researchers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)