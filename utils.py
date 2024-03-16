import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import os
from collections import Counter
from Bio.SeqUtils import nt_search
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
rna_fold_path = os.path.expanduser('~/.conda/envs/myenv/bin/RNAfold')
os.environ["PATH"] += os.pathsep + os.path.dirname(rna_fold_path)
import json
import RNA
import itertools
from tqdm import tqdm
import pandas as pd
# ----------------- General utility functions -----------------

def read_fasta(file_path):
    """Read a FASTA file and return a dictionary with sequence IDs as keys and sequences as values."""
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, 'fasta')}

def rnafold(seq):
    """Run RNAfold to get the structure of a sequence."""
    result = subprocess.run(['RNAfold'], input=seq, text=True, capture_output=True)
    # Parse the RNAfold output to get just the structure string
    structure = result.stdout.split('\n')[1].split(' ')[0]
    return structure

def get_ct_data(sequence, structure, output_file):
    ct_data = {}
    with open(output_file, 'w') as f:
        f.write(f"{len(sequence)} ENERGY\n")
        stack = []
        
        for i, (nt, struct) in enumerate(zip(sequence, structure), start=1):
            prev_nt = i - 1 if i != 1 else 0
            next_nt = i + 1 if i != len(sequence) else 0
            
            if struct == '(':
                stack.append(i)
                ct_data[i] = [i, nt, prev_nt, next_nt, 0, i]  # Temporarily set as unpaired
            elif struct == ')':
                partner = stack.pop()
                f.write(f"{i} {nt} {prev_nt} {next_nt} {partner} {i}\n")
                ct_data[i] = [i, nt, prev_nt, next_nt, partner, i]
                ct_data[partner][4] = i  # Update partner's pairing information
            else:
                f.write(f"{i} {nt} {prev_nt} {next_nt} 0 {i}\n")
                ct_data[i] = [i, nt, prev_nt, next_nt, 0, i]

    return ct_data

def load_specific_columns(csv_path):
    """
    Load a CSV file and return a DataFrame with specific columns.

    Parameters:
    - csv_path (str): The path to the CSV file.

    Returns:
    - pd.DataFrame: A DataFrame containing only the specified columns.
    """
    # Attempt to read the CSV file
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file at {csv_path} was not found.")
    except Exception as e:
        raise Exception(f"An error occurred while reading {csv_path}: {e}")

    # Specify the required columns
    required_columns = ["full_seq", "full_seq_folding", "Mature", "Star"]

    # Check if all required columns are in the DataFrame
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in the CSV file: {missing_columns}")
        
    # Filter out rows where any of the required columns contain non-string (specifically float, as it likely indicates NaNs)
    for col in required_columns:
        df = df[pd.to_numeric(df[col], errors='coerce').isna()]  # Keeps rows where conversion to numeric fails (i.e., strings)

    df[required_columns] = df[required_columns].astype(str)

    # Remove rows with 'nan' values if they were converted to string 'nan'
    df = df[~df[required_columns].isin(['nan']).any(axis=1)]

    # Return the DataFrame with only the required columns
    return df[required_columns]

def dotbracket_to_ct(sequence, structure):
    # Initialize the list for CT data
    ct_data = []
    
    # Get the pair table from the dot-bracket structure
    # The pair table is 1-indexed and the first element is the length of the RNA
    pair_table = RNA.ptable(structure)
    
    # Loop through each nucleotide in the sequence
    for i, nt in enumerate(sequence, start=1):
        # Get the pairing partner from the pair table
        partner = pair_table[i]
        
        # Write the CT data: index, nt, prev, next, partner, index
        prev_nt = i - 1 if i > 1 else 0
        next_nt = i + 1 if i < len(sequence) else 0
        
        # Append the data as a tuple or a list
        ct_data.append((i, nt, prev_nt, next_nt, partner, i))
    
    return ct_data

                


# ----------------- mirGenDb Data Creation ----------------- #

def prepare_mirgendb_data(full_seq_file, mature_file, star_file):
    # Read the sequences from the FASTA files
    full_seqs = read_fasta(full_seq_file)
    matures = read_fasta(mature_file)
    stars = read_fasta(star_file)
    
    # Prepare the DataFrame
    df = pd.DataFrame.from_dict(full_seqs, orient='index', columns=['full_seq'])
    
    # Add mature and star sequences by matching IDs
    df['Mature'] = df.index.map(matures).fillna('')
    df['Star'] = df.index.map(stars).fillna('')
    tqdm.pandas()
    # Compute RNA folding for full sequences
    
    print('adding rna fold - take long time')
    df['full_seq_folding'] = df['full_seq'].progress_apply(rnafold)
    
    return df

# ----------------- mirMachine Data Creation ----------------- #

def build_dict():
    """
    Build and return an empty dictionary with predefined keys.
    
    :return: Empty dictionary with keys
    """
    keys = ['Chr', 'Start_hairpin', 'End_hairpin', 'Strand', 'Hairpin_seq', 'Mature_connections', 'Mature_BP_ratio',
            'Mature_max_bulge', 'Loop_length', 'Valid mir', 'Loop_seq', 'Mature', 'Mature_Length', 'Mature2', '3p/5p',
            'Hairpin_seq_trimmed', 'Star', 'Start_star', 'End_star', 'Star_length', 'Star_connections', 'Star_BP_ratio',
            'Star_max_bulge', 'Hairpin_seq_trimmed_length', 'Window', 'Max_bulge_symmetry', 'full_seq',
            'full_seq_folding',
            'seed', 'flank1', 'flank2', 'UG', 'UGUG']
    return {i: [] for i in keys}

def parse_gff_file(gff_file_path, seed_family_dict):
    with open(gff_file_path) as file:
        lines = file.readlines()

    # Extract metadata
    metadata = {}
    for line in lines:
        if line.startswith("# ") and not line.startswith("##"):
            try:
                key, value = line.strip("# ").split(":", 1)
                metadata[key.strip()] = value.strip()
            except ValueError:
                continue

    genome_file = metadata.get('Genome file', 'Unknown')
    species = metadata.get('Species', 'Unknown')
    
    data = []
    filtered_count = 0 ; more_than_once = 0
    for line in lines:
        if line.startswith("#") or line.startswith("\n"):
            continue
        try:
            parts = line.split("\t")
            attributes = dict(attr.split("=") for attr in parts[-1].split(";") if "=" in attr)
            gene_id = attributes.get("gene_id", "").replace(".PRE", "")
            full_seq = attributes.get("sequence_with_30nt", "").replace('\n','')
            hairpin_trimmed = full_seq[30:-30]
            # print(hairpin_trimmed
            seeds = [seed for seed in seed_family_dict.get(gene_id.upper(), []) if seed in hairpin_trimmed]
            if seeds:
                chosen_seed = seeds[0]  # Select the first seed if multiple are found
                seed_start = hairpin_trimmed.find(chosen_seed)
                seed_end = seed_start + len(chosen_seed)
                mature = hairpin_trimmed[max(0, seed_start-1):min(len(hairpin_trimmed), seed_end+14)]
                if hairpin_trimmed.count(chosen_seed) == 1:
                    data.append([full_seq, gene_id, genome_file, species, hairpin_trimmed, chosen_seed, mature])  # Adding count as 1
                else:
                    more_than_once +=1
            else:
                filtered_count +=1
            # If no seed is found, do not add the row to data or count in hairpin greater than 1 
        except Exception as e:
            print(f"Error processing line: {line}\nError: {e}")
    if data:
        df = pd.DataFrame(data, columns=['full_seq', 'gene_id', 'Genome file', 'Species', 'hairpin_trimmed', 'seed', 'mature'])
    else:
        # Return an empty DataFrame with the expected columns if no data was added
        df = pd.DataFrame(columns=['full_seq', 'gene_id', 'Genome file', 'Species', 'hairpin_trimmed', 'seed', 'mature'])

    # No need to filter out rows without a seed as they were never added to the data list
    return df, (filtered_count,more_than_once)

def run_features_on_sebastian_db(gff_path, output_path, seed_family_csv_path):
    # Read seed family dictionary
    features_dict = build_dict()
    seed_family_df = pd.read_csv(seed_family_csv_path, encoding='ISO-8859-1')
    seed_family_df['Seed'] = seed_family_df['Seed'].apply(lambda x: x.replace('U', 'T') if isinstance(x, str) else x)
    seed_family_dict = seed_family_df.groupby('Family')['Seed'].apply(list).to_dict()

    all_dfs = []
    more_than_once = 0 ; filtered_count =0
    for i, gff_filename in enumerate(os.listdir(gff_path)):
        if gff_filename.endswith('.gff'):
            gff_file_path = os.path.join(gff_path, gff_filename)
            df, counts = parse_gff_file(gff_file_path, seed_family_dict)
            more_than_once += counts[1] ;filtered_count+=counts[0]
            all_dfs.append(df)
            
            if i % 100 == 0:
                print(f'{i} gffs processed.')
    
    # Combine all dataframes
    gff_data = pd.concat(all_dfs, ignore_index=True)
    gff_data.to_csv(output_path, index=False)
    print(f"Processing complete. Output saved to {output_path}")
    return (more_than_once,filtered_count)

def filter_mirMachine_data(mirMachine_path, mirGendb_path): # First filter of mirMachine data

    # Load the CSV files
    bastian_df = pd.read_csv(mirMachine_path)
    mirgene_db = pd.read_csv(mirGendb_path)

    # Add length columns
    bastian_df['mature_len'] = bastian_df['mature'].apply(lambda x: len(x))
    bastian_df['hairpin_trimmed_len'] = bastian_df['hairpin_trimmed'].apply(lambda x: len(x))
    bastian_df['full_seq_len'] = bastian_df['full_seq'].apply(lambda x: len(x))

    # Analyze the `hairpin_trimmed` column lengths and filter
    hairpin_lengths = bastian_df['hairpin_trimmed'].str.len()
    filtered_by_length = bastian_df[(hairpin_lengths >= 45) & (hairpin_lengths <= 65)]

    # Filter based on sequence presence in `miRGeneDB`
    mirgene_sequences = set(mirgene_db['full_seq'])
    filtered_by_presence = filtered_by_length[~filtered_by_length['full_seq'].isin(mirgene_sequences)]

    # Filter based on mature sequence length
    filtered_by_mature_length = filtered_by_presence[filtered_by_presence['mature'].str.len() == 22]

    # Count rows filtered by each criterion
    count_length_filtered = len(bastian_df) - len(filtered_by_length)
    count_presence_filtered = len(filtered_by_length) - len(filtered_by_presence)
    count_mature_length_filtered = len(filtered_by_presence) - len(filtered_by_mature_length)

    # Save filtered counts and reasons to a text file
    with open('./filtered_counts.txt', 'w') as f:
        f.write(f"Rows filtered by length: (less than 45 or more than 65): {count_length_filtered}\n")
        f.write(f"Rows filtered by presence in miRGeneDB: {count_presence_filtered}\n")
        f.write(f"Rows filtered by mature sequence not being 22: {count_mature_length_filtered}\n")
    print("Filtering complete. Results saved to 'filtered_counts.txt'")
    # print summery of lengths
    # len(filtered_by_length) , len(filtered_by_presence) , len(bastian_df) - len(filtered_by_length), len(filtered_by_mature_length)
    print(f"Length of filtered_by_length: {len(filtered_by_length)}")
    print(f"Length of filtered_by_presence: {len(filtered_by_presence)}")
    print(f"Length of bastian_df - filtered_by_length: {len(bastian_df) - len(filtered_by_length)}")
    print(f"Length of filtered_by_mature_length: {len(filtered_by_mature_length)}")
    return filtered_by_mature_length


def find_star_sequence(hairpin, mature):
    # Compute the secondary structure
    ss, _ = RNA.fold(hairpin)
    pairing_list = RNA.ptable(ss)
    mature_start = hairpin.find(mature)

    # Determine the direction of the mature sequence
    if mature_start <= 5:
        direction = '5p'
    else:
        direction = '3p'

    # Initialize variables
    star_start = None
    star_end = None

    # Logic for finding the STAR sequence start and end
    if direction == '5p':
        # Find the end of the STAR sequence for mature at the start
        paired_pos = pairing_list[mature_start + len(mature)]
        if paired_pos != 0:
            # STAR sequence starts two nucleotides after the paired position
            star_start = paired_pos + 2 -1 
        else:
            # Find the next paired nucleotide if the last of mature is not paired
            for i in range(len(mature) - 1, -1, -1):
                if pairing_list[mature_start + i] != 0:
                    dist = len(mature) - i 
                    star_start = pairing_list[mature_start + i] + 2 - dist - 1
                    break
    else:
        # Mature is at the end, find the start of the STAR sequence
        paired_pos = pairing_list[mature_start]
        if paired_pos != 0:
            # STAR sequence ends two nucleotides before the paired position
            star_end = paired_pos + 2 -1 
        else:
            # Find the next paired nucleotide if the first of mature is not paired
            for i in range(len(mature)):
                if pairing_list[mature_start + i] != 0:
                    star_end = pairing_list[mature_start + i] - 2 + i -1 
                    break

    # Adjusting the logic for extracting the STAR sequence based on direction and the new logic
    star_sequence = None
    if direction == '5p' and star_start is not None:
        star_sequence = hairpin[star_start:]
    elif direction == '3p' and star_end is not None:
        star_sequence = hairpin[:star_end]
        star_start = 0  # STAR starts at the beginning of the hairpin for 3p direction

    if star_sequence:
        star_length = len(star_sequence)
    else:
        return None

    return star_sequence, star_start, star_length, ss, mature_start

def add_star_and_folding_information(df):
    """
    Applies additional processing to a DataFrame by adding star sequence, star start position,
    star sequence length, hairpin folding structure, and mature start position for each row.

    Parameters:
    - df: DataFrame to process.
    
    find_star_sequence: A function that takes hairpin_trimmed and mature sequence as inputs
                          and returns the star sequence, its start position, length, the hairpin folding structure,
                          and mature start position.

    Returns:
    - DataFrame: The input DataFrame with additional columns for star sequence information and hairpin folding.
    """

    from tqdm.auto import tqdm

    # Initialize new columns with None values
    df['star'] = None
    df['star_start'] = None
    df['star_length'] = None
    df['hairpin_folding'] = None
    df['mature_start'] = None

    # Process each row in the DataFrame
    for index, row in tqdm(df.iterrows(), total=df.shape[0]):
        results = find_star_sequence(row['hairpin_trimmed'], row['mature'])
        if results:
            star_sequence, star_start, star_length, fold, mature_start = results
            df.at[index, 'star'] = star_sequence
            df.at[index, 'star_start'] = int(star_start) + 1
            df.at[index, 'star_length'] = star_length
            df.at[index, 'hairpin_folding'] = fold
            df.at[index, 'mature_start'] = int(mature_start) + 1

    return df

def filter_dataset(df, star_length_threshold=(22,25), mature_length_threshold=(22,23), loop_length_threshold=(8,25)):
    # Filter based on the given thresholds
    filtered_df = df[
        (df['star_length'] >= star_length_threshold[0]) & (df['star_length'] <= star_length_threshold[1]) &
        (df['mature_len'] >= mature_length_threshold[0]) & (df['mature_len'] <= mature_length_threshold[1]) &
        (df['loop_length'] >= loop_length_threshold[0]) & (df['loop_length'] <= loop_length_threshold[1])
    ]

    return filtered_df


def update_mature_sequences(row):
    if pd.notna(row['mature_start']) and row['mature_start'] > 30:
        new_mature = row['hairpin_trimmed'][row['mature_start']-1:]
        row['mature'] = new_mature
        row['mature_len'] = len(new_mature)
    return row


def calculate_loop_length(row):
    # Assuming the mature sequence comes after the star sequence
    # This will need adjustment if the assumption doesn't always hold
    if row['star_start'] < row['mature_start']:  # Star comes first
        loop_start = row['star_start'] + row['star_length']
        loop_end = row['mature_start'] - 1
    else:  # Mature comes first
        loop_start = row['mature_start'] + row['mature_len']
        loop_end = row['star_start'] - 1
    return max(0, loop_end - loop_start + 1)

def final_mirMachine_preprocess(df,save_path = None):

    print("Running final preprocessing on mirMachine data")
    df = add_star_and_folding_information(df)
    print("Star sequence and folding information added")
    
    df = df.apply(update_mature_sequences, axis=1)
    print("Mature sequences updated")
    
    df['loop_length'] = df.apply(calculate_loop_length, axis=1)

    filtered_df = filter_dataset(df)
    print("Dataset filtered based on lengths thresholds")
    
    filtered_df['full_seq_folding'] = [run_rnafold(x) for x in tqdm(filtered_df['full_seq'].tolist(), desc="Running RNAfold")]
    print("Secondary structure calculated for all mirMachine sequences")
    
    filtered_df = filtered_df.rename(columns={"mature": "Mature", "star": "Star"})
    
    filtered_df_with_features = extract_features_from_df(filtered_df)
    print("Finish extracting features from the mirMachine data")

    # Save the filtered DataFrame with star and other updated information
    if save_path:
        filtered_df_with_features.to_csv(save_path, index=False)
        print(f"Final Mirmachine data saved to {save_path}")
    else:
        return filtered_df_with_features
    
    return filtered_df_with_features


# ----------------- Exracting Features ----------------- #

def calculate_statistics(sequences):
    num_sequences = len(sequences)
    avg_length = int(sum(len(seq) for seq in sequences) / num_sequences)
    min_length = min(len(seq) for seq in sequences)
    max_length = max(len(seq) for seq in sequences)

    return {
        "Number of sequences": num_sequences,
        "Average sequence length": avg_length,
        "Min sequence length": min_length,
        "Max sequence length": max_length
    }
    
# Define Translation Dictionary
translation_dict = {
    '(C': 'P', ')C': 'Q', '.C': 'E',
    '(G': 'R', ')G': 'H', '.G': 'I',
    '(T': 'J', ')T': 'K', '.T': 'L',
    '(A': 'M', ')A': 'N', '.A': 'O',
}

def encode_sequence(row):
    # Get Full Sequence and Folding Sequence from the row
    full_seq = row['full_seq']
    full_seq_folding = row['full_seq_folding']
    
    # Encode based on the Translation Dictionary
    encoded_seq = ''.join([translation_dict.get(full_seq_folding[i] + full_seq[i], full_seq[i]) for i in range(len(full_seq))])
    
    # Encode Star and Mature sequences with [ZZZZZ,BBBB,DDDDD,FFFFF] special tokens
    mature_start = full_seq.find(row['Mature'])
    mature_end = mature_start + len(row['Mature'])
    star_start = row['Start_star']
    star_end = row['End_star']
    
    if star_start < star_end:
        if mature_start < star_start:
            encoded_seq = (encoded_seq[:mature_start] + 'ZZZZZ' +
                        encoded_seq[mature_start:mature_end] + 'BBBBB' +
                        encoded_seq[mature_end:star_start] + 'DDDDD' +
                        encoded_seq[star_start:star_end] + 'FFFFF' +
                        encoded_seq[star_end:])
        else:
            encoded_seq = (encoded_seq[:star_start] + 'DDDDD' +
                        encoded_seq[star_start:star_end] + 'FFFFF' +
                        encoded_seq[star_end:mature_start] + 'ZZZZZ' +
                        encoded_seq[mature_start:mature_end] + 'BBBBB' +
                        encoded_seq[mature_end:])
    
    return encoded_seq

    # Define Reverse Translation Dictionary
reverse_translation_dict = {
    'P': 'C', 'Q': 'C', 'E': 'C',
    'R': 'G', 'H': 'G', 'I': 'G',
    'J': 'T', 'K': 'T', 'L': 'T',
    'M': 'A', 'N': 'A', 'O': 'A',
}


def calculate_global_alignment(sequence1, sequence2):
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    return alignments[0][2]


def calculate_novelty_global(sequences, training_set):
    total_score = 0
    for generated in tqdm(sequences):
        # Calculate the average global alignment score of a generated sequence to all sequences in the training set
        scores = [calculate_global_alignment(generated, train_seq) for train_seq in training_set]
        avg_score = sum(scores) / len(scores)
        total_score += avg_score
    return total_score / len(sequences)

def calculate_local_diversity(gen_seq, real_sequences):
    local_scores = [(pairwise2.align.localxx(gen_seq, real_seq, score_only=True), real_seq) for real_seq in real_sequences]
    max_score, max_real_seq = max(local_scores, key=lambda x: x[0])
    normalized_max_score = max_score / min(len(gen_seq), len(max_real_seq))
    avg_score = sum(score for score, _ in local_scores) / len(real_sequences)
    return (round(avg_score, 2), round(normalized_max_score, 2), max_real_seq)

def calculate_global_diversity(gen_seq, real_sequences):
    global_scores = [(pairwise2.align.globalxx(gen_seq, real_seq, score_only=True), real_seq) for real_seq in real_sequences]
    max_score, max_real_seq = max(global_scores, key=lambda x: x[0])
    normalized_max_score = max_score / min(len(gen_seq), len(max_real_seq))
    avg_score = sum(score for score, _ in global_scores) / len(real_sequences)
    return (round(avg_score, 2), round(normalized_max_score, 2), max_real_seq)

# Helper function to count connections (in mature or star)
def count_connections(start, end, ct_data):
    connections = 0
    for i in range(start, end + 1):
        if ct_data[i][4] != 0 and not (start <= ct_data[i][4] <= end):
            connections += 1
    return connections

def calculate_max_bulge(mature_range, star_range, ct_data):
    max_bulge_mature = 0
    max_bulge_star = 0
    current_bulge = 0
    in_bulge = False
    bulge_start = -1
    bulge_end = -1
    mature_bulges = []
    star_bulges = []

    for i in range(1, len(ct_data)):
        if ct_data[i][4] == 0:  # If current position is unpaired
            if not in_bulge:
                bulge_start = i-1  # Start of the bulge
                # print(i,bulge_start)
                in_bulge = True
            current_bulge += 1
            # print(ct_data[i+1] , ct_data[i+1][4])
            if i == len(ct_data)-1 or ct_data[i+1][4] != 0:  # If it's the last position or next is paired
                bulge_end = i+1  # End of the bulge
                # print(i,bulge_end)
                # Check if the bulge is completely contained within the mature range
                if bulge_start >= mature_range[0] and bulge_end <= mature_range[1]:
                    # print(f"Mature Bulge found from index {bulge_start} to {bulge_end}. Size: {current_bulge}")
                    mature_bulges.append((bulge_start, bulge_end))
                    max_bulge_mature = max(max_bulge_mature, current_bulge)
                # Check if the bulge is completely contained within the star range
                elif bulge_start >= star_range[0] and bulge_end <= star_range[1]:
                    # print(f"Star Bulge found from index {bulge_start} to {bulge_end}. Size: {current_bulge}")
                    star_bulges.append((bulge_start, bulge_end))
                    max_bulge_star = max(max_bulge_star, current_bulge)
                current_bulge = 0
                in_bulge = False
        else:
            current_bulge = 0
            in_bulge = False

    return {
        "mature_max_bulge": max_bulge_mature,
        "star_max_bulge": max_bulge_star,
        "mature_bulges": mature_bulges,
        "star_bulges": star_bulges
    }


def debug_print(debug, *args):
    if debug:
        print(*args)

def check_seed_family(seed):
    # Load seed families from CSV

    filename = os.path.abspath("/sise/vaksler-group/IsanaRNA/Transformers/Rom/Data_source/seed_family_from_mirgendb.csv")

    # Read CSV using pandas
    df = pd.read_csv(filename, encoding='ISO-8859-1')

    # Create a dictionary from the 'Seed' and 'Family' columns
    seed_family_dict = df.set_index('Seed')['Family'].to_dict()
    
    return seed_family_dict.get(seed, "Unknown")

# Helper function to calculate UG & UGUG
def find_ug_sequences(decoded_seq, mature_start, mature_end, threshold=3):
    # Calculate the search ranges considering the threshold
    ug_search_start = max(mature_end - 14 - threshold, 0)
    ug_search_end = min(mature_start + threshold, len(decoded_seq))
    
    ugug_search_start = max(mature_end + 1 - threshold, 0)
    ugug_search_end = min(mature_end + 3 + threshold, len(decoded_seq))
    
    # Search for 'UG' and 'UGUG' sequences in the calculated ranges
    ug_index = decoded_seq.find('UG', ug_search_start, ug_search_end) + 1  # +1 to make it 1-indexed
    ug = ug_index if ug_index != 0 else "FALSE"
    
    ugug_index = decoded_seq.find('UGUG', ugug_search_start, ugug_search_end) + 1  # +1 to make it 1-indexed
    ugug = ugug_index if ugug_index != 0 else "FALSE"
    
    return ug, ugug

# Helper function to mer features (for sizes 1 or 2)
def calculate_mer_ratios(sequence, mer_size):
    total_length = len(sequence)
    # If mer_size is greater than the sequence length, return an error or handle appropriately
    if mer_size > total_length:
        return "(0.00 .. 0.00)"
    
    # Generate all possible combinations of nucleotides of length mer_size
    possible_mers = [''.join(p) for p in itertools.product('AUCG', repeat=mer_size)]
    mer_counts = Counter([sequence[i:i+mer_size] for i in range(total_length - mer_size + 1)])
    
    # Calculate ratios
    mer_ratios = {mer: mer_counts[mer] / total_length for mer in possible_mers}
    max_ratio = max(mer_ratios.values(), default=0)
    min_ratio = min(mer_ratios.values(), default=0)
    
    return f"({min_ratio:.2f} .. {max_ratio:.2f})"

def calculate_energy(sequence):
    # Calculate the secondary structure and the free energy of the structure
    structure, energy = RNA.fold(sequence)
    
    return energy

def is_valid_dot_bracket(structure, sequence):
    stack = []
    valid_pairs = [('U', 'G'), ('G', 'U'), ('C', 'G'), ('G', 'C'), ('A', 'U'), ('U', 'A')]
    invalid_pairs_dict = {}
    
    for i, char in enumerate(structure):
        if char == '(':
            stack.append((i, sequence[i]))
        elif char == ')':
            if not stack:
                return False, invalid_pairs_dict
            else:
                position, nucleotide = stack.pop()
                if (nucleotide, sequence[i]) not in valid_pairs:
                    invalid_pair = f"{nucleotide}{sequence[i]}"
                    invalid_pairs_dict[invalid_pair] = invalid_pairs_dict.get(invalid_pair, 0) + 1
                    
    # Check if there are any unmatched parentheses left
    if stack:
        return False, invalid_pairs_dict
    
    # Check if any invalid pairs were found
    if invalid_pairs_dict:
        return False, invalid_pairs_dict
        
    return True, invalid_pairs_dict

def run_rnafold(rna_sequence):
    # Start the RNAfold process
    process = subprocess.Popen(['RNAfold'],
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               text=True)

    # Send the RNA sequence and get the output
    output, error = process.communicate(rna_sequence)

    # Check for errors
    if error:
        print("Error:", error)
    else:
        return output

def retrieve_parts_only_nts(encoded_seq):

    # Identify positions of special tokens
    mature_start = encoded_seq.find('ZZZZZ')
    end_mature_indices = encoded_seq[mature_start + 5:].find('BBBBB')
    mature_end = mature_start + 5 + end_mature_indices

    star_start = encoded_seq.find('DDDDD')
    end_star_indices = encoded_seq[star_start + 5:].find('FFFFF')
    star_end = star_start + 5 + end_star_indices
    
    # Decode mature and star parts
    decoded_mature = encoded_seq[mature_start + 5:mature_end].replace('T', 'U')
    decoded_star = encoded_seq[star_start + 5:star_end].replace('T', 'U')
    
    decoded_mature = decoded_mature.replace('ZZZZZ', '').replace('BBBBB', '').replace('DDDDD', '').replace('FFFFF', '')
    decoded_star = decoded_star.replace('ZZZZZ', '').replace('BBBBB', '').replace('DDDDD', '').replace('FFFFF', '')
    
    intermediate_seq = encoded_seq.replace('ZZZZZ', '').replace('BBBBB', '').replace('DDDDD', '').replace('FFFFF', '')

    # Check for -1 in the find results
    if mature_start == -1 or end_mature_indices == -1 or star_start == -1 or end_star_indices == -1:
        return None

    # Decode full sequence and folding
    full_seq_folding = run_rnafold(intermediate_seq)
    if full_seq_folding:
        # print(full_seq_folding) 
        full_seq_folding = full_seq_folding.split('\n')[-2][:-9] # extract nesscery structre from rnafold output
        # print(len(full_seq_folding), len(intermediate_seq)) 
    else:
        return -1
    full_seq = intermediate_seq.replace('T', 'U')


    return {
        'encoded_seq': encoded_seq,
        'decoded_seq': full_seq,
        'mature': decoded_mature,
        'star': decoded_star,
        'full_seq_folding': full_seq_folding
    }


def adjust_index(index):
    return index + 1

def extract_features_only_nts(decoded_seq, full_seq_folding, mature, star):
    # Save CT file
    if is_valid_dot_bracket(full_seq_folding,decoded_seq):
        ct_data = get_ct_data(decoded_seq, full_seq_folding, "output.ct")
    else:
        return None

    # Locate mature and star sequences within the full sequence
    mature_start = decoded_seq.find(mature)
    mature_end = mature_start + len(mature) - 1  # Adjust to 0-based indexing
    star_start = decoded_seq.find(star)
    star_end = star_start + len(star) - 1  # Adjust to 0-based indexing
    if mature_start < 0 or star_start < 0:
        return None

    if mature_start < star_start:
        loop = decoded_seq[mature_end+1:star_start]
        flank1 = decoded_seq[:mature_start]
        flank2 = decoded_seq[star_end+1:]
        direction = '5p'
    else:
        loop = decoded_seq[star_end+1:mature_start]
        flank1 = decoded_seq[:star_start]
        flank2 = decoded_seq[mature_end+1:]
        direction = '3p'
    # Compute the seed section
    seed_start = mature_start + 2  # 2nd nucleotide, 0-indexed
    seed_end = seed_start + 7  # 7th nucleotide, 0-indexed
    seed = decoded_seq[seed_start-1:seed_end-1]  
    seed_family = check_seed_family(seed) # loading dict from shared folder 
    

    # Calculate the indices for Loop, flank1, and flank2
    loop_start_index = decoded_seq.find(loop) + 1
    loop_end_index = loop_start_index + len(loop) - 1

    flank1_start_index = decoded_seq.find(flank1) + 1
    flank1_end_index = flank1_start_index + len(flank1) - 1

    flank2_start_index = decoded_seq.find(flank2) + 1
    flank2_end_index = flank2_start_index + len(flank2) - 1
    
    
    mature_range, star_range = (adjust_index(mature_start), adjust_index(mature_end)), (adjust_index(star_start), adjust_index(star_end))
    
    # Feature calculations
    # print(mature_start, mature_end , len(ct_data))
    mature_connections = count_connections(adjust_index(mature_start), adjust_index(mature_end), ct_data)
    mature_bp_ratio = mature_connections / len(mature) if len(mature) != 0 else 0
    calculate_max_bulge_result = calculate_max_bulge(mature_range, star_range, ct_data)
    mature_max_bulge, star_max_bulge = calculate_max_bulge_result["mature_max_bulge"], calculate_max_bulge_result["star_max_bulge"]
    mature_bulges = calculate_max_bulge_result["mature_bulges"]
    star_bulges = calculate_max_bulge_result["star_bulges"]
    # mature_max_asymmetry = calculate_bulge_asymmetry(mature_bulges, ct_data , 'mature')
    # print('mature',mature_bulges)
    # star_max_asymmetry = calculate_bulge_asymmetry(star_bulges, ct_data,'star')
    # print('star',star_bulges)
    star_connections = count_connections(adjust_index(star_start), adjust_index(star_end), ct_data)
    star_bp_ratio = star_connections / len(star) if len(star) != 0 else 0
    ug, ugug = find_ug_sequences(decoded_seq, mature_start, mature_end, threshold=10)
    hairpin_trimmed = decoded_seq.replace(flank1, '').replace(flank2,'')
    h_start = adjust_index(decoded_seq.find(hairpin_trimmed))
    h_end = mature_start + len(hairpin_trimmed)
    energy = calculate_energy(hairpin_trimmed)
    one_mer_mature = calculate_mer_ratios(mature, 1)
    two_mer_mature = calculate_mer_ratios(mature, 2)
    one_mer_h_trimm = calculate_mer_ratios(hairpin_trimmed, 1)
    two_mer_h_trimm = calculate_mer_ratios(hairpin_trimmed, 2)
    one_mer_full = calculate_mer_ratios(decoded_seq, 1)
    two_mer_full = calculate_mer_ratios(decoded_seq, 2)
    # ct_data2 = dotbracket_to_ct(decoded_seq, full_seq_folding)


    feature_dict = {
        'full_seq': decoded_seq,
        'Mature': mature,
        'Mature_Length': len(mature),
        'Star': star,
        'End_star': adjust_index(star_end),
        'full_seq_folding': full_seq_folding,
        'Loop_seq': loop,
        'Loop_length': len(loop),
        'flank1': flank1,
        'flank2': flank2,
        'Mature_start': adjust_index(mature_start),
        'Mature_end': adjust_index(mature_end),
        'Star_start': adjust_index(star_start),
        'Star_end': adjust_index(star_end),
        'Star_length': len(star),
        'Mature_length': len(mature),
        'Loop_seq_start': loop_start_index,
        'Loop_seq_end': loop_end_index,
        'flank1_start': flank1_start_index,
        'flank1_end': flank1_end_index,
        'flank2_start': flank2_start_index,
        'flank2_end': flank2_end_index,
        'seed': seed,
        'seed_start': seed_start,
        'seed_end': seed_end,
        'seed_family': seed_family,
        '3p/5p': direction,
        'Mature_connections': mature_connections,
        'Mature_BP_ratio': round(mature_bp_ratio, 2),
        'Mature_max_bulge': mature_max_bulge,
        'Star_connections': star_connections,
        'Star_BP_ratio': round(star_bp_ratio, 2),
        'Star_max_bulge': star_max_bulge,
        'UG': ug,
        'UGUG': ugug,
        'hairpin_trimmed': hairpin_trimmed,
        'hairpin_trimmed_length': len(hairpin_trimmed),
        'folding_energy': round(energy, 2),
        'one_mer_mature': one_mer_mature,
        'two_mer_mature': two_mer_mature,
        'one_mer_hairpin_trimmed': one_mer_h_trimm ,
        'two_mer_hairpin_trimmed': two_mer_h_trimm ,
        'one_mer_full': one_mer_full,
        'two_mer_full': two_mer_full,

    }

    return feature_dict

def process_rna_sequences(rna_sequences):
    results = []
    invalid_results = []
    
    summary = {
        'multiple_mature': 0,
        'multiple_star': 0,
        'missing_mature': 0,
        'missing_star': 0,
        'invalid_dot_bracket': 0,
        'invalid_pairs': {},
        'total_sequences': len(rna_sequences),
        'unique_before_filtering': len(set(rna_sequences)),  # Number of unique sequences before filtering
        'unique_after_filtering': 0,  # To be calculated after filtering
        'unique_invalid': 0,
        'failed_rnafold': 0  
    }
    
    for seq in tqdm(rna_sequences, desc="Processing sequences"):
        if 'ZZZZZ' not in seq or 'BBBBB' not in seq:
            summary['missing_mature'] += 1
        elif seq.count('ZZZZZ') > 1 or seq.count('BBBBB') > 1:
            summary['multiple_mature'] += 1

        if 'DDDDD' not in seq or 'FFFFF' not in seq:
            summary['missing_star'] += 1
        elif seq.count('DDDDD') > 1 or seq.count('FFFFF') > 1:
            summary['multiple_star'] += 1

        parts = retrieve_parts_only_nts(seq)  # Ensure retrieve_parts is defined somewhere
        if parts:
            if parts == -1:
                summary['failed_rnafold'] += 1
                continue
            is_valid, invalid_pairs = is_valid_dot_bracket(parts['full_seq_folding'], parts['decoded_seq'])
            if is_valid:
                results.append(parts)
            else:
                summary['invalid_dot_bracket'] += 1
                invalid_result = {'sequence': parts['full_seq_folding'], 'encoded_seq': seq,'invalid_pairs':invalid_pairs}
                invalid_result.update(invalid_pairs)
                invalid_results.append(invalid_result)
                for pair, count in invalid_pairs.items():
                    summary['invalid_pairs'][pair] = summary['invalid_pairs'].get(pair, 0) + count

    results_df = pd.DataFrame(results).drop_duplicates(subset='decoded_seq')
    invalid_results_df = pd.DataFrame(invalid_results)
    
    # Calculate unique values after filtering
    summary['unique_after'] = results_df['decoded_seq'].nunique() if 'decoded_seq' in results_df.columns else 0
    summary['unique_invalid'] = invalid_results_df['sequence'].nunique() if 'sequence' in invalid_results_df.columns else 0
    
    return results_df, invalid_results_df, summary

def extract_features_from_df(rna_df):
    features = []
    # Add tqdm around your loop
    print('Starting extract features')
    for index, row in tqdm(rna_df.iterrows(), total=len(rna_df), desc="Extracting features"):
        feature = extract_features_only_nts(row['full_seq'], row['full_seq_folding'], row['Mature'], row['Star'])
        if feature:
            features.append(feature)

    features_df = pd.DataFrame(features)
    return features_df

def add_diversity_to_df(rna_df, real_sequences):
    # Calculate local and global diversity scores
    tqdm.pandas() # This is for using tqdm with pandas' apply method

    rna_df['avg_local_diversity'], rna_df['max_local_diversity'], rna_df['max_real_seq_local_alignment'] = zip(*rna_df['full_seq'].progress_apply(lambda x: calculate_local_diversity(x, real_sequences)))

    rna_df['avg_global_diversity'], rna_df['max_global_diversity'], rna_df['max_real_seq_global_alignment'] = zip(*rna_df['full_seq'].progress_apply(lambda x: calculate_global_diversity(x, real_sequences)))

    return rna_df


