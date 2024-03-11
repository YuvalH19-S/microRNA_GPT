# microRNA_GPT

This repository contains scripts and tools for collecting, filtering, and extracting biological features from microRNA data, as well as for training GPT models on this data.

## Data Collection and Preprocessing

### Preprocessing Pipeline Overview

The `preprocess.py` script contain the data collection, filtering, and feature extraction process. It gets data paths configurations specified in `preprocess_config.ini` for easy modification of input and output paths without needing to alter the script directly.

### `utils.py`

The `utils.py` file contains essential functions used by `preprocess.py`, including data collection routines, filtering mechanisms, and feature extraction algorithms.

### Configuration with `preprocess_config.ini`

Paths for input data and where to save the output data

### Running the Preprocessing Pipeline

1. Configure `preprocess_config.ini` with the appropriate paths and settings for your data and processing needs.
2. Run the pipeline:
    ```bash
    python preprocess.py
    ```
    This will collect data as specified, apply filters, extract features, and output the processed datasets ready for machine learning model training.

### Tokenization
After that, we run in **colab** the code that is in the file 'tokenization.ipynb', where we create and train the tokenizer on the data (bastian + mirGeneDB) to use later.

### Training notebooks
Then, we use 'GPT_pretrained_mature_star_after_preprocess.ipynb', this script use the tokenizer we've trained, (with or without flanks), then we train a model based on bastian data and then second train on mirGeneDB data.

After that, we can split the data to train and test. We do this using 'preprocess_cluster.ipynb', where we remove duplicates of pre-mirna and create clusters, based on mature similarity > 80%. We take some of the clusters to test - we use inly the human data from this cluseters to be the test.

The script we can use now is 'Human_fine_tune_star_mature.ipynb', there we do fine tune on the previous model, using human data (with or without flanks). We can generate sequences - full seq or with completions of mature/star - using the test data from the previous script.

Lastly, we use the script 'Extract_features_only_nts.ipynb' to extract features from the generated sequences.

## [Placeholder for Training Scripts]

Training scripts and methodologies will be added here. These scripts will detail the process for training GPT models on the prepared microRNA dataset, including model configuration, training parameters, and evaluation techniques.

## Getting Started

(Instructions on how to use the scripts, set up the environment, etc.)

## Contributing

(Guidelines for contributing to the repository.)


