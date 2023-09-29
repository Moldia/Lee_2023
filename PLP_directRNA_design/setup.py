from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

VERSION = '0.1.1'
DESCRIPTION = 'Package to design Padlock probes based on a reference transcriptome'
LONG_DESCRIPTION = 'Package to design Padlock probes based on a reference transcriptome'

# Setting up
setup(
    name="PLP_directRNA_design",
    version=VERSION,
    author="Sergio Marco Salas",
    author_email="<sergiomarco.salas@scilifelab.se>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=['pandas>=1.1.5','Bio>=1.0','numpy>=1.19.0'],
    keywords=['python', 'spatial transcriptomics', 
            'spatial resolved transcriptomics', 
            'in situ sequencing', 
            'ISS','padlock probe design'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Researchers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
