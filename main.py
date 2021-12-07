import pandas as pd
from scripts.ProteinComplexSimulation import ProteinComplexes
from scripts.CompareComplexes import CompareProteinComplexes
from scripts.Filters import Filter
from scripts.ReporterVisualizer import ReporterVisualizer
import yaml
from yaml import CLoader as Loader
from scripts.utils import create_folder
import logging


def main(input1,input2,settings={}):
    
    create_folder(settings)

    if len(settings['FILTER SETTINGS']['FILTER'])>0:
        logging.info("Filtering data...")
        filt1=Filter(settings)
        input1=filt1.run_filters(input1)
        input2=filt1.run_filters(input2)
        
    logging.info("Comparing data...")
    cc=CompareProteinComplexes(input1,input2,settings)
    cc.run_comparisons()
    logging.info("Creating reports and data visualizations...")
    rv=ReporterVisualizer(compared_complex=cc,settings=settings)
    rv.run_all()


if __name__ == "__main__":
    settings=None
    with open('./settings.yaml','r') as f:
            settings=yaml.load(f,Loader=Loader)

    input1=pd.read_csv(f'./data/{settings["INPUT1"]["FILENAME"]}',delimiter=settings["INPUT1"]["DELIMITER"])
    input2=pd.read_csv(f'./data/{settings["INPUT2"]["FILENAME"]}',delimiter=settings["INPUT2"]["DELIMITER"])

    main(input1,input2,settings=settings)
