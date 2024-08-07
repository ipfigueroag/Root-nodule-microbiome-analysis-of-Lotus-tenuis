import os
import re

def extract_substring(line):
    # This regex pattern captures the part of the string after 'rawdata.' and before '_1.fq.gz'
    match = re.search(r'rawdata\.([A-Za-z0-9_]+)_\d\.fq\.gz', line)
    if match:
        return match.group(1)
    else:
        return None


def process_file(input_file, output_file):
    unique_extracted = set() # set to track unique substrings
    data_dict = {} # Dictionary to store the corresponding lines for _1 and _2 terminations

    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip() # Remove leading and trailing whitespace
            extracted = extract_substring(line)
            if extracted:
                if extracted not in data_dict:
                    data_dict[extracted] = {}
                absolute_path = os.path.abspath(line) # Get the absolute path of the file
                if line.endswith('_1.fq.gz'):
                    data_dict[extracted]['_1'] = absolute_path
                elif line.endswith('_2.fq.gz'):
                    data_dict[extracted]['_2'] = absolute_path
    
    with open(output_file, 'w') as outfile:
        for key in data_dict:
            if key not in unique_extracted:
                if '_1' in data_dict[key] and '_2' in data_dict[key]:
                    unique_extracted.add(key)
                    outfile.write(f"{key}\t{data_dict[key]['_1']}\t{data_dict[key]['_2']}\n")

input_file = 'filenames.txt'  # Replace with your input file path
output_file = 'manifest.txt'  # Replace with your desired output file path

process_file(input_file, output_file)

