mkdir ISS_bioinfo
cd ISS_bioinfo

git clone -b dev_branch https://github.com/Moldia/ISS_preprocessing
cd ISS_preprocessing
conda env create --file 'preprocessing.yml' --name 1_ISS_preprocessing
conda activate 1_ISS_preprocessing
python setup.py install
pip install ipykernel
python -m ipykernel install --user --name="1_ISS_preprocessing"
conda deactivate
cd ..

git clone https://github.com/Moldia/ISS_decoding
cd ISS_decoding
conda env create --file 'decoding.yml' --name 2_ISS_decoding
conda activate 2_ISS_decoding
python setup.py install
python -m ipykernel install --user --name= "2_ISS_decoding"
conda deactivate
cd ..

git clone https://github.com/Moldia/ISS_postprocessing
cd ISS_postprocessing
conda env create --file 'postprocessing.yml' --name 3_ISS_postprocessing
conda activate 3_ISS_postprocessing
python setup.py install
pip install stardist
pip install ipykernel
python -m ipykernel install --user --name= "3_ISS_postprocessing"
conda deactivate
cd ..

