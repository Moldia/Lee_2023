import pciSeq
from pciSeq import utils
import pandas as pd


def preprocess_spots(spots_file, conversion_factor = 0.1625):
    spots = pd.read_csv(spots_file)
    spots = iss_spots.dropna()
    spots = spots[['target', 'xc', 'yc']]
    spots = spots.rename(columns = {'target':'Gene', 'xc': 'x', 'yc': 'y'})
    spots['x'] = spots['x']/conversion_factor
    spots['y'] = spots['y']/conversion_factor
    return spots

def get_most_probable_call_pciseq(cellData):
    # create dataframe with most probable cell types
    append_data = pd.DataFrame(columns=['Cell_Num', 'X', 'Y', 'ClassName', 'Prob'])
    # create dataframe with most probable cell types
    for i, cell in enumerate(cellData.Cell_Num):
        cell_names = cellData.ClassName[i]
        cell_prob = cellData.Prob[i]
        try:
            max_prob = max(cell_prob)
        except TypeError: 
            max_prob = max([cell_prob])
        try: 
            index = [i for i, j in enumerate(cell_prob) if j == max_prob]
        except TypeError: 
            index = [i for i, j in enumerate([cell_prob]) if j == max_prob]
        cellname = cell_names[index[0]]
        X = cellData.X[i]
        Y = cellData.Y[i]
        data = [cell, X, Y, cellname, max_prob]
        append_data.loc[i] = data
    return append_data

def run_pciseq(spots, coo_mask, sc_expression_matrix, output_dir, save_output = True):
    cellData, geneData = pciSeq.fit(spots, coo_mask, sc_expression_matrix)
    cellData = cellData.reset_index()
    most_probable = get_most_probable_call_pciseq(cellData)
    if save_output == True:
        cellData.to_json(output_dir+'/cellData.json')
        geneData.to_json(output_dir+'/geneData.json')
        most_probable.to_csv(output_dir+'/most_probable.csv')
    return cellData, geneData, most_probable