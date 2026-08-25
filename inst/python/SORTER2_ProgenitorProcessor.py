from Bio import SeqIO
import pandas as pd
import argparse
import os
from collections import defaultdict

from sorter2_progenitor_match import match_progenitor_row

parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-indir", "--inputdir", default=None)
parser.add_argument("-map", "--mapfile")
parser.add_argument("-og", "--outgroups")
parser.add_argument("-min", "--minseq")
parser.add_argument("-dif", "--filterundiff")
parser.add_argument("-v", "--verbose", action="store_true", default=False)
args = parser.parse_args()

# -indir defaults to -wd for backward compatibility
indir = args.inputdir if args.inputdir else args.workingdir

print('Collapsing progenitors to clade sets')

def rename_fasta_sequence_ids(input_folder, output_folder, csv_file):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_file)

    # Ensure the CSV file has the correct columns
    if df.shape[1] != 3:
        raise ValueError(
            "CSV file must have exactly three columns: [hybrid, progenitor, clade]"
        )

    # Convert the CSV data into a list of dictionaries for easy lookup
    csv_data = df.to_dict(orient='records')

    # Species-level (hybrid = bare species code) mapfile rows that
    # actually got used, collected so the monophyly warning prints
    # once per species rather than once per matched record.
    species_level_matches = set()

    # Process each FASTA file in the input folder
    for fasta_filename in os.listdir(input_folder):
        if fasta_filename.endswith(".fasta"):
            fasta_path = os.path.join(input_folder, fasta_filename)
            new_fasta_path = os.path.join(
                output_folder,
                fasta_filename.replace(".fasta", "_collapsed.fasta")
            )

            with open(new_fasta_path, "w") as output_handle:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    seq_id_parts = record.id.split('_')

                    clade, warn_species_key = match_progenitor_row(seq_id_parts, csv_data)
                    if clade is not None:
                        seq_id_parts[2] = clade
                    if warn_species_key is not None:
                        species_level_matches.add(warn_species_key)

                    # Create the new sequence ID
                    new_seq_id = '_'.join(seq_id_parts)
                    record.id = new_seq_id
                    record.description = new_seq_id

                    # Write the modified record to the new FASTA file
                    SeqIO.write(record, output_handle, "fasta")

            if args.verbose:
                print(f"Processed {fasta_filename} -> {new_fasta_path}")

    for species_key in sorted(species_level_matches):
        print(
            f"Warning: mapfile row(s) for hybrid={species_key!r} matched by "
            f"species code alone and applied to every voucher with this "
            f"code. If {species_key!r} is not monophyletic, specify "
            f"individual vouchers instead (e.g. "
            f"'<voucher>_{species_key}')."
        )

# -map and -og are now full paths
mapfile = args.mapfile

rename_fasta_sequence_ids(indir, args.workingdir, mapfile)

print('Removing sequences that mapped to outgroups')

def remove_outgroups(fasta_folder, text_file):
    # Read the list of strings to match from the text file
    with open(text_file, 'r') as file:
        remove_list = set(line.strip() for line in file.readlines())

    # Process each FASTA file in the folder
    for fasta_filename in os.listdir(fasta_folder):
        if fasta_filename.endswith("_collapsed.fasta"):
            fasta_path = os.path.join(fasta_folder, fasta_filename)
            new_fasta_path = os.path.join(
                fasta_folder,
                fasta_filename.replace("_collapsed.fasta", "_cladecollapse1.fasta")
            )

            with open(new_fasta_path, "w") as output_handle:
                records_to_keep = []
                for record in SeqIO.parse(fasta_path, "fasta"):
                    seq_id_parts = record.id.split('_')

                    # Check if the third underscore-delimited item is in the remove list
                    if seq_id_parts[2] not in remove_list:
                        records_to_keep.append(record)

                # Write the filtered records to the new FASTA file
                SeqIO.write(records_to_keep, output_handle, "fasta")

            if args.verbose:
                print(f"Processed {fasta_filename} -> {new_fasta_path}")

# Remove outgroup sequences
outgroupfile = args.outgroups
remove_outgroups(args.workingdir, outgroupfile)

# Remove intermediate files
for fasta in os.listdir(args.workingdir):
    if fasta.endswith("_collapsed.fasta"):
        os.remove(os.path.join(args.workingdir, fasta))

def filter_min_sequences(fasta_folder, min_count):
    # Initialize a dictionary to track the occurrence of third items by unique sample ID
    sample_third_item_counts = defaultdict(lambda: defaultdict(int))
    sequences_to_remove = defaultdict(set)

    # First pass: count occurrences of third underscore-delimited items
    for fasta_filename in os.listdir(fasta_folder):
        if fasta_filename.endswith("_cladecollapse1.fasta"):
            fasta_path = os.path.join(fasta_folder, fasta_filename)

            for record in SeqIO.parse(fasta_path, "fasta"):
                seq_id_parts = record.id.split('_')

                # Ensure the sequence ID has at least four underscore-delimited items
                if len(seq_id_parts) >= 4:
                    sample_id = '_'.join(seq_id_parts[:2])
                    third_item = seq_id_parts[2]

                    # Increment the count of the third item for this sample ID
                    sample_third_item_counts[sample_id][third_item] += 1

    # Identify sequences to remove
    for sample_id, third_items in sample_third_item_counts.items():
        for third_item, count in third_items.items():
            if count < min_count:
                sequences_to_remove[sample_id].add(third_item)

    # Second pass: remove sequences with third items below the threshold
    for fasta_filename in os.listdir(fasta_folder):
        if fasta_filename.endswith("_cladecollapse1.fasta"):
            fasta_path = os.path.join(fasta_folder, fasta_filename)
            new_fasta_path = os.path.join(
                fasta_folder,
                fasta_filename.replace("_cladecollapse1.fasta", "_cladecollapse.fasta")
            )

            with open(new_fasta_path, "w") as output_handle:
                records_to_keep = []

                for record in SeqIO.parse(fasta_path, "fasta"):
                    seq_id_parts = record.id.split('_')

                    # Ensure the sequence ID has at least four underscore-delimited items
                    if len(seq_id_parts) >= 4:
                        sample_id = '_'.join(seq_id_parts[:2])
                        third_item = seq_id_parts[2]

                        # Keep the record only if its third item is not in the removal list
                        if third_item not in sequences_to_remove[sample_id]:
                            records_to_keep.append(record)
                    else:
                        # If the sequence ID doesn't have four parts, keep it as is
                        records_to_keep.append(record)

                # Write the filtered records to the new FASTA file
                SeqIO.write(records_to_keep, output_handle, "fasta")

            if args.verbose:
                print(f"Processed {fasta_filename} -> {new_fasta_path}")

# Filter Sequences
minseqfilter = int(int(args.minseq)*2)
filter_min_sequences(args.workingdir, minseqfilter)

# Remove intermediate files
for fasta in os.listdir(args.workingdir):
    if fasta.endswith("collapse1.fasta"):
        os.remove(os.path.join(args.workingdir, fasta))

# Remove duplicate sequences

print('Keeping longest sequence if two identical sequence ids are present ''\n')

def remove_dup(input_file, output_file):
    sequences = list(SeqIO.parse(input_file, "fasta"))
    sequence_dict = {}

    for seq in sequences:
        seq_id = seq.id
        if seq_id not in sequence_dict:
            sequence_dict[seq_id] = seq
        else:
            # Compare sequences with the same ID and keep the longest one with the most base pairs
            if len(seq.seq) > len(sequence_dict[seq_id].seq):
                sequence_dict[seq_id] = seq

    # Write the filtered sequences to the output FASTA file
    SeqIO.write(sequence_dict.values(), output_file, "fasta")

for file in os.listdir(args.workingdir):
    if file.endswith('_cladecollapse.fasta'):
        fp = os.path.join(args.workingdir, file)
        remove_dup(fp, fp[:-6] + '_duprem.fasta')
        os.remove(fp)

if args.filterundiff == 'T':

    # Define function to filter sequences that don't have atleast 2 unique progenitor names represented in locus sequences

    def filter_sequences(input_file):
        # Read the input FASTA file
        records = list(SeqIO.parse(input_file, "fasta"))

        # Group sequences based on the first two tab-delimited items
        grouped_sequences = {}
        for record in records:
            header_items = record.description.split("_")
            sample_id = tuple(header_items[:2])

            if len(header_items) == 4:
                # Get sequences with 4 tab-delimited items
                if sample_id not in grouped_sequences:
                    grouped_sequences[sample_id] = set()

                grouped_sequences[sample_id].add(header_items[2])

        # Identify hybrid to be removed that only has one progenitor represented at the locus
        samples_to_remove = {
            sample_id
            for sample_id, unique_strings in grouped_sequences.items()
            if len(unique_strings) < 2
        }

        # Filter sequences
        filtered_records = [
            record for record in records
            if len(record.description.split("_")) == 3
            or tuple(record.description.split("_")[:2]) not in samples_to_remove
        ]

        # Write the filtered sequences to a new FASTA file
        output_file = input_file.replace("_duprem.fasta", "_differentiated.fasta")
        SeqIO.write(filtered_records, output_file, "fasta")


    for file in os.listdir(args.workingdir):
        if file.endswith('cladecollapse_duprem.fasta'):
            fp = os.path.join(args.workingdir, file)
            filter_sequences(fp)
            os.remove(fp)

# Generate summary table of hybrid progenitor distributions across loci

def analyze_progenitor_distributions(fasta_files, sample_names, output_csv):
    try:
        # Extract samples identifiers from sample names
        with open(sample_names, 'r') as file:
            sample_ids_list = [line.strip() for line in file]

        sample_ids = {sample_id: None for sample_id in sample_ids_list}
        if args.verbose:
            print("Sample IDs extracted: ", sample_ids)
        # Make dictionary to hold the unique count data
        data = {os.path.basename(fasta): {sample: set() for sample in sample_ids.keys()} for fasta in fasta_files}
        if args.verbose:
            print("Initialized data structure: ", data)

        for fasta_file in fasta_files:
            fasta_name = os.path.basename(fasta_file)
            if args.verbose:
                print(f"Processing file: {fasta_file}")
            # Read the FASTA file
            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = record.id
                # Extract the first two underscore-delimited strings
                sample_id = '_'.join(seq_id.split('_')[:2])

                if sample_id in sample_ids:
                    # Extract the third underscore-delimited string
                    unique_str = seq_id.split('_')[2]

                    if fasta_name not in data:
                        data[fasta_name] = {}

                    if sample_id not in data[fasta_name]:
                        data[fasta_name][sample_id] = set()

                    data[fasta_name][sample_id].add(unique_str)

        # Convert sets to counts
        for fasta_name in data:
            for sample in data[fasta_name]:
                data[fasta_name][sample] = len(data[fasta_name][sample])

        # Create a DataFrame from the data dictionary
        df = pd.DataFrame(data).T
        df.index.name = 'alignment'
        df.reset_index(inplace=True)

        # Write the DataFrame to a CSV file
        df.to_csv(output_csv, index=False)

    except Exception as e:
        print(f"An error occurred: {e}")

# Get Hybrid sample IDs
def find_unique_sample_ids(fasta_folder, output_file):
    unique_sample_ids = set()  # Use a set to store unique sample IDs

    # Iterate through each FASTA file in the specified folder
    for fasta_filename in os.listdir(fasta_folder):
        if "cladecollapse_" in fasta_filename:
            fasta_path = os.path.join(fasta_folder, fasta_filename)

            # Parse the FASTA file
            for record in SeqIO.parse(fasta_path, "fasta"):
                seq_id_parts = record.id.split('_')

                # Check if the sequence ID has at least 4 underscore-delimited items
                if len(seq_id_parts) >= 4:
                    # Create the unique sample ID by joining the first two parts with an underscore
                    sample_id = '_'.join(seq_id_parts[:2])
                    # Add the unique sample ID to the set
                    unique_sample_ids.add(sample_id)

    # Write the unique sample IDs to the output file
    with open(output_file, 'w') as f:
        for sample_id in sorted(unique_sample_ids):
            f.write(sample_id + '\n')

    if args.verbose:
        print(f"Found {len(unique_sample_ids)} unique sample IDs. Output written to {output_file}")

# Get Hybrid IDs
find_unique_sample_ids(args.workingdir, os.path.join(args.workingdir, 'hybrid_ids.txt'))

# Get FASTA file names
fasta_files = [
    os.path.join(args.workingdir, file)
    for file in os.listdir(args.workingdir)
    if "cladecollapse_" in file
]

sample_names = os.path.join(args.workingdir, 'hybrid_ids.txt')

output_csv = os.path.join(args.workingdir, "progenitor_distributions.csv")

analyze_progenitor_distributions(fasta_files, sample_names, output_csv)

# Summarize progenitor clades across all loci

def summarize_progenitors(fasta_files, sample_names, output_csv):
    try:
        with open(sample_names, 'r') as file:
            sample_ids_list = [line.strip() for line in file]

        sample_ids = {sample_id: None for sample_id in sample_ids_list}
        if args.verbose:
            print("Sample IDs extracted: ", sample_ids)
        # Make a dictionary to hold the unique count data
        data = {sample: [] for sample in sample_ids.keys()}

        for fasta_file in fasta_files:
            fasta_name = os.path.basename(fasta_file)
            if args.verbose:
                print(f"Processing file: {fasta_file}")
            # Temporary dictionary for this FASTA file
            temp_data = {sample: set() for sample in sample_ids.keys()}

            # Read the FASTA file
            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_id = record.id
                # Extract the first two underscore-delimited strings
                sample_id = '_'.join(seq_id.split('_')[:2])

                if sample_id in sample_ids:
                    # Extract the third underscore-delimited string
                    unique_str = seq_id.split('_')[2]
                    temp_data[sample_id].add(unique_str)

            # Convert sets to counts and store in main data dictionary
            for sample in temp_data:
                data[sample].append(len(temp_data[sample]))

        # Determine the maximum number of unique strings for column headers
        all_counts = [c for counts in data.values() for c in counts]
        if not all_counts:
            print("No hybrid samples detected; skipping progenitor summary.")
            return
        max_unique_count = max(all_counts)

        # Initialize a dictionary for the final counts
        final_counts = {
            sample: {f'{i}_progenitors': 0 for i in range(max_unique_count + 1)}
            for sample in sample_ids.keys()
        }

        # Populate the final counts dictionary
        for sample, counts in data.items():
            for count in counts:
                final_counts[sample][f'{count}_progenitors'] += 1

        # Create a DataFrame from the final counts dictionary
        df = pd.DataFrame(final_counts).T
        df.index.name = 'samples'
        df.reset_index(inplace=True)

        # Write the DataFrame to a CSV file
        df.to_csv(output_csv, index=False)

    except Exception as e:
        print(f"An error occurred: {e}")

summary_csv = os.path.join(args.workingdir, "progenitor_summary.csv")

summarize_progenitors(fasta_files, sample_names, summary_csv)
