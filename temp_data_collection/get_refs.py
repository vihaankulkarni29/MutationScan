import urllib.request
import os

# UniProt Accession IDs for E. coli K-12 MG1655 targets
targets = {
    "acrA_WT.faa": "P0AE06",
    "acrB_WT.faa": "P31224",
    "TolC_WT.faa": "P02930",
    "acrR_WT.faa": "P0AEL6",
    "marR_WT.faa": "P0A878"
}

os.makedirs("data/refs", exist_ok=True)

for filename, uniprot_id in targets.items():
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    filepath = os.path.join("data/refs", filename)
    print(f"Downloading {filename}...")
    try:
        urllib.request.urlretrieve(url, filepath)
    except Exception as e:
        print(f"Failed to download {filename}: {e}")

print("All Wild-Type references downloaded successfully!")
