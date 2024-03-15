# preprocess.py
import configparser
import logging
from utils import run_features_on_sebastian_db,load_specific_columns,extract_features_from_df,filter_mirMachine_data,final_mirMachine_preprocess

def main():
    # Setup logging
    logging.basicConfig(filename='preprocessing.log', level=logging.INFO, 
                        format='%(asctime)s:%(levelname)s:%(message)s')

    # Load configurations
    config = configparser.ConfigParser()
    config.read('preprocess_config.ini')

    # Paths
    seed_family_csv = str(config['Paths']['seed_family_csv'])
    mirMachine_path = str(config['Paths']['mirMachine_path'])
    mirGendb_path = str(config['Paths']['mirGendb_path'])

    # Constant Paths
    mirMachine_DB = config['Constant Paths']['mirMachine_DB']
    star_mirgeneDB = config['Constant Paths']['star_mirgeneDB']
    mature_mirgeneDB = config['Constant Paths']['mature_mirgeneDB']
    fullseq_mirgeneDB = config['Constant Paths']['fullseq_mirgeneDB']
    mirgen_with_features = config['Constant Paths']['mirgeneDB_star_mature']

    # Processing Steps
    try:
        logging.info("Loading raw data from MirMachine DB")
        counts = run_features_on_sebastian_db(mirMachine_DB, mirMachine_path, seed_family_csv)
        logging.info(f"Total mirMachine from raw data: {counts}")

        logging.info("Preparing MirGeneDB data -  Takes long time")
        mirgen_db = load_specific_columns(mirgen_with_features)
        mirgen_db_features = extract_features_from_df(mirgen_db)
        mirgen_db_features.to_csv(mirGendb_path, index=False)
        logging.info("MirGeneDB updated with features")

        logging.info("First filter on mirMachine data")
        mirmachine_after_filters = filter_mirMachine_data(mirMachine_path, mirGendb_path)

        logging.info("Starting final preprocess on MirMachine data")
        final_mirMachine = final_mirMachine_preprocess(mirmachine_after_filters, save_path=mirMachine_path)
        logging.info("Final preprocess completed")

    except Exception as e:
        logging.error("Error during preprocessing: " + str(e))

    logging.info("Preprocessing pipeline completed successfully")

if __name__ == "__main__":
    main()
