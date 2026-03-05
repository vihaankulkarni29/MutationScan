import pandas as pd
from pathlib import Path


def clean_and_filter_date(date_str):
    if pd.isna(date_str) or date_str == 'N/A' or date_str == 'missing':
        return 'Unknown'
    try:
        year = int(str(date_str)[:4])
        if year == 1905:  # Rescue known THSTI 2023 typo
            return 'Rescued (Typo)'
        if year >= 2015:
            return 'Pass'
        return 'Too Old'
    except Exception:
        return 'Unknown'


def main():
    root_input = Path('pure_indian_metadata.csv')
    temp_input = Path('temp_data_collection/pure_indian_metadata.csv')
    input_file = root_input if root_input.exists() else temp_input

    df = pd.read_csv(input_file)
    df['Date_Check'] = df['Collection_Date'].apply(clean_and_filter_date)
    final_df = df[df['Date_Check'].isin(['Pass', 'Rescued (Typo)', 'Unknown'])]

    print(f"Original Strains: {len(df)}")
    print(f"Ancient Strains Dropped: {len(df[df['Date_Check'] == 'Too Old'])}")
    print(f"Final Modern Dataset: {len(final_df)}")

    final_df.drop(columns=['Date_Check']).to_csv('final_clinical_dataset.csv', index=False)


if __name__ == '__main__':
    main()