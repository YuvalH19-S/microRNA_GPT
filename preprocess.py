# preprocess.py
import configparser
import random
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import subprocess
import itertools
from collections import Counter
import logging
from utils import *

def main():
    # Setup logging
    logging.basicConfig(filename='preprocessing.log', level=logging.INFO, 
                        format='%(asctime)s:%(levelname)s:%(message)s')

    # Load configurations
    config = configparser.ConfigParser()
    config.read('preprocess_config.ini')

    # Paths
    seed_family_csv = config['Paths']['seed_family_csv']
    mirMachine_path = config['Paths']['mirMachine_path']
    mirGendb_path = config['Paths']['mirGendb_path']

    # Constant Paths
    mirMachine_DB = config['Constant Paths']['mirMachine_DB']
    star_mirgeneDB = config['Constant Paths']['star_mirgeneDB']
    mature_mirgeneDB = config['Constant Paths']['mature_mirgeneDB']
    fullseq_mirgeneDB = config['Constant Paths']['fullseq_mirgeneDB']

    # Processing Steps
    try:
        logging.info("Loading raw data from MirMachine DB")
        counts = run_features_on_sebastian_db(mirMachine_DB, mirMachine_path, seed_family_csv)
        logging.info(f"Total mirMachine from raw data: {counts}")

        logging.info("Preparing MirGeneDB data")
        mirgen_db = prepare_mirgendb_data(fullseq_mirgeneDB, mature_mirgeneDB, star_mirgeneDB)
        mirgen_db_features = extract_features_from_df(mirgen_db)
        mirgen_db_features.to_csv(mirGendb_path, index=False)
        logging.info("MirGeneDB features saved")

        logging.info("First filter on mirMachine data")
        mirmachine_after_filters = filter_mirMachine_data(mirMachine_path, mirgen_db_features)

        logging.info("Starting final preprocess on MirMachine data")
        final_mirMachine = final_mirMachine_preprocess(mirmachine_after_filters, save_path=mirMachine_path)
        logging.info("Final preprocess completed")

    except Exception as e:
        logging.error("Error during preprocessing: " + str(e))

    logging.info("Preprocessing pipeline completed successfully")

if __name__ == "__main__":
    main()
