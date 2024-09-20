import pandas as pd
import numpy as np
import sys
from scipy.stats import gmean

def replace_underscores(df):
    # Iterate over each column in the DataFrame
    for col in df.columns:
        # Check if the column's data type is object, typically used for strings
        if df[col].dtype == 'object':
            df[col] = df[col].str.replace('_', ' ', regex=True)
    return df
    
def extract_needed_numeric_data(file, group_ids, ignore_ids):
    # Load csv
    df = pd.read_csv(file)

    # Split the group_ids string into a list of columns
    group_columns = group_ids.split(',')
    ignore_columns = ignore_ids.split(',')
     
    # Remove columns not needed (ignore_ids)
    df = df.drop(columns=ignore_columns, errors='ignore')

    # Replace underscores for latex
    df = replace_underscores(df)

    # Select only numeric columns, but ensure group columns are included
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()  # Get list of numeric column names
    all_relevant_cols = list(set(numeric_cols + group_columns))  # Combine and remove duplicates

    # Filter DataFrame to include only the relevant columns (numeric + group columns)
    df_relevant = df[all_relevant_cols]

    return df_relevant, group_columns
    
def compute_mean(file, group_ids, ignore_ids, outfile):

    df_numeric, group_columns = extract_needed_numeric_data(file, group_ids, ignore_ids)
    
    # Compute mean of each numeric column grouped by the specified group columns
    df_grouped = df_numeric.groupby(group_columns).agg(
        lambda x: round(np.mean(x.dropna()), 2))

    df_grouped.to_csv(outfile)
    

def compute_gmean(file, group_ids, ignore_ids, outfile):

    df_numeric, group_columns = extract_needed_numeric_data(file, group_ids, ignore_ids)

    # Compute geometric mean of each numeric column grouped by the specified group columns
    df_grouped = df_numeric.groupby(group_columns).agg(lambda x: round(gmean(x.dropna() if x.count() > 0 else np.nan),2))

    df_grouped.to_csv(outfile)

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 6:
        print("Usage: python script.py <file> <group_id> <ignore_id> <arithmetic(a)/geometric(g)> <outfilename>")
        sys.exit(1)

    # Unpack command line arguments
    file = sys.argv[1]
    group_id = sys.argv[2]
    ignore_id = sys.argv[3]
    which_mean = sys.argv[4]
    outfile = sys.argv[5]

    # Run the aggregation function
    if which_mean == "a":
    	compute_mean(file, group_id, ignore_id, outfile+"_mean.csv")
    else:
        if which_mean == "g":
	        compute_gmean(file, group_id, ignore_id, outfile+"_gmean.csv")
        else:
            print("specify mean: a for arithmetic mean and g for geometric mean")
