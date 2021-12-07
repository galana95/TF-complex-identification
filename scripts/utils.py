import logging
import sys
import os
from datetime import datetime
def is_equal(string1,string2):
    if isinstance(string1,list):string1=string1[0]
    if isinstance(string2,list):string2=string2[0]
    equal=False
    if string1.lower() == string2.lower(): equal=True
    return equal
def logprint(string,level='info'):
    if level == 'info':
        logging.info(string)
    elif level == 'warning':
        logging.warning(string)
    elif level == 'critical':
        logging.critical(string)
def get_folder_name(settings):
    dtime=datetime.today().strftime('%Y%m%d')
    i1_name=settings['INPUT1']['NAME']
    i2_name=settings['INPUT2']['NAME']
    i1i2_name=i1_name[:5].lower()+'_'+i2_name[:5].lower()
    run_name=settings['RUN NAME']
    directory=f'./results/{run_name}_{i1i2_name}_{dtime}'
    return directory
def create_folder(settings):
    directory=get_folder_name(settings)


    overwrite=settings['REPORTER VISUALIZER SETTINGS']['OVERWRITE EXISTING']
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.makedirs(directory+"/complex_target_subnetworks")
        os.makedirs(directory+"/figures")
        logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(f'{directory}/logging.log'),
            logging.StreamHandler()
        ])
    else:
        logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(f'{directory}/logging.log'),
            logging.StreamHandler()
        ])
        logging.warning(f'{directory} already exists. Overwriting is set to {overwrite}. Abort if necessary.')
    