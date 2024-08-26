import pandas as pd

# This script was created with assistance from the AI tool ChatGPT and was manually refined and checked for quality.

# Adjusted function to analyze datasets and return results per replicate
def analyze_dataset(file_paths, total_reads_list):
    all_replicates_data = []

    for file_path, total_reads in zip(file_paths, total_reads_list):
        conversion_counts = {}
        total_coverage = 0

        with open(file_path, 'r') as file:
            next(file)  # Skip the header line
            for line in file:
                columns = line.strip().split('\t')
                coverage = int(columns[4])
                base_count = eval(columns[6])  # Convert string representation of list to an actual list

                total_coverage += coverage

                for i, count in enumerate(base_count):
                    if count > 0:
                        for j, other_count in enumerate(base_count):
                            if i != j and other_count > 0:
                                conversion_pair = f"{['A', 'C', 'G', 'T'][i]}->{['A', 'C', 'G', 'T'][j]}"
                                conversion_counts[conversion_pair] = conversion_counts.get(conversion_pair, 0) + ((other_count * coverage) / total_reads)

        # Convert counts to decimal format normalized by total genome coverage
        normalized_conversion_rates = {pair: count / total_coverage for pair, count in conversion_counts.items()}
        all_replicates_data.append((file_path, normalized_conversion_rates))

    return all_replicates_data

# Write results to a CSV file
def write_results_to_csv(data, filename):
    with open(filename, 'w') as csv_file:
        csv_file.write("File,Conversion Pair,Rate\n")
        for file_path, rates in data:
            for pair, rate in rates.items():
                csv_file.write(f"{file_path},{pair},{rate:.6f}\n")

# Analyze datasets and print results
def print_and_save_results(file_paths, total_reads_list, label):
    results = analyze_dataset(file_paths, total_reads_list)
    print(f"\nNormalized Conversion Rates ({label}):")
    for file_path, rates in results:
        print(f"\n{file_path}:")
        for pair, rate in rates.items():
            print(f"{pair}: {rate:.6f}")
    write_results_to_csv(results, f"{label}_normalized_conversion_rates.csv")

# Running the analysis for both datasets
print_and_save_results(file_paths_control, total_reads_control, "Control")
print_and_save_results(file_paths_experimental, total_reads_experimental, "Experimental")
