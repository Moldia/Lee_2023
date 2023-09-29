from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.0.23'
DESCRIPTION = 'Package used to decode preprocessed ISS images, including SpaceTX formatting and plotting'
LONG_DESCRIPTION = 'This package is used to decode ISS data from appropriately preprocessed images.'

# Setting up
setup(
    name="ISS_decoding",
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