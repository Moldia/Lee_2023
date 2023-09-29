from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.0'
DESCRIPTION = 'Package used to preprocess and process ISS data'
LONG_DESCRIPTION = 'This package can be used to process ISS data, that includes the preprocessing of image data (aligning images across squencing cycles and the subsequent stitching of tiles and tiling). In addition, there are functions included which allows for the formatting of the data to SpaceTX format and the subsequent decoding.'

# Setting up
setup(
    name="ISS_processing",
    version=VERSION,
    author="Christoffer Mattsson Langseth",
    author_email="<christoffer.langseth@scilifelab.se>",
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
