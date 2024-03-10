# microRNA_GPT

This repository contains scripts and tools for collecting, filtering, and extracting biological features from microRNA data, as well as for training GPT models on this data.

## Data Collection and Preprocessing

The initial script, `preapare_data_gff.ipynb`, is designed to collect microRNA data from various sources, apply a series of filters to ensure data quality, and extract relevant biological features for further analysis. This script lays the foundation for generating a high-quality dataset suitable for training machine learning models, including GPT-based architectures, to predict or analyze microRNA sequences.

After that, we run in **colab** the code that is in the file 'tokenization.ipynb', where we create and train the tokenizer on the data (bastian + mirGeneDB) to use later.

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


