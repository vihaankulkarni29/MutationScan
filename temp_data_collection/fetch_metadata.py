import pandas as pd
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import time
import json

# Replace with your actual email and API key
EMAIL = "vihaankulkarni29@gmail.com"
API_KEY = "c734d921ce3764ef09fcab1c1499e4f41508"

def esearch_biosample(query):
    """Searches the BioSample database and returns the ID."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'biosample',
        'term': f'"{query}"[Strain]',
        'retmode': 'json',
        'email': EMAIL,
        'api_key': API_KEY
    }
    url = f"{base_url}?{urllib.parse.urlencode(params)}"
    try:
        response = urllib.request.urlopen(url)
        data = json.loads(response.read().decode('utf-8'))
        id_list = data['esearchresult'].get('idlist', [])
        return id_list[0] if id_list else None
    except Exception as e:
        print(f"Error searching for {query}: {e}")
        return None

def efetch_biosample_details(biosample_id):
    """Hunts the XML for exact city, hospital, and GPS coordinates."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'biosample',
        'id': biosample_id,
        'retmode': 'xml',
        'email': EMAIL,
        'api_key': API_KEY
    }
    url = f"{base_url}?{urllib.parse.urlencode(params)}"
    
    # Initialize with N/A
    metadata = {
        'Exact_Location': 'N/A', 
        'GPS_Lat_Lon': 'N/A',
        'Hospital_or_Institute': 'N/A',
        'Collection_Date': 'N/A', 
        'Host': 'N/A', 
        'Isolation_Source': 'N/A',
        'Clinical_Disease': 'N/A'
    }
    
    try:
        response = urllib.request.urlopen(url)
        xml_data = response.read()
        root = ET.fromstring(xml_data)
        
        # Aggressively parse all Attributes
        for attribute in root.findall('.//Attribute'):
            attr_name = attribute.get('attribute_name', '').lower()
            text_val = attribute.text if attribute.text else 'N/A'
            
            if 'geo_loc_name' in attr_name:
                metadata['Exact_Location'] = text_val
            elif 'lat_lon' in attr_name:
                metadata['GPS_Lat_Lon'] = text_val
            elif 'collected_by' in attr_name or 'isolate_name_alias' in attr_name:
                metadata['Hospital_or_Institute'] = text_val
            elif 'collection_date' in attr_name:
                metadata['Collection_Date'] = text_val
            elif 'host' in attr_name:
                metadata['Host'] = text_val
            elif 'isolation_source' in attr_name:
                metadata['Isolation_Source'] = text_val
            elif 'host_disease' in attr_name or 'disease' in attr_name:
                metadata['Clinical_Disease'] = text_val
                
        return metadata
    except Exception as e:
        print(f"Error fetching details for ID {biosample_id}: {e}")
        return metadata

def main():
    print("Loading BVBRC genome data...")
    df = pd.read_csv('BVBRC_genome_amr.csv')
    
    strains = df['Genome Name'].dropna().unique()
    print(f"Found {len(strains)} unique strains. Preparing Full Geospatial Interrogation...\n")
    
    results = []
    
    # Process ALL strains (removed the [:10] slice)
    for i, strain in enumerate(strains):
        print(f"[{i+1}/{len(strains)}] Tracking Strain: {strain}")
        
        clean_strain = strain.replace("Escherichia coli strain ", "").replace("Escherichia coli ", "").strip()
        
        bs_id = esearch_biosample(clean_strain)
        if bs_id:
            meta = efetch_biosample_details(bs_id)
            print(f"  └─ Location: {meta['Exact_Location']}")
            
            row = {'Strain': strain, 'BioSample_ID': bs_id}
            row.update(meta)
            results.append(row)
        else:
            print("  └─ Not found in primary NCBI BioSample index.")
            
        time.sleep(0.15) # Safe API pacing

    # Create master dataframe
    output_df = pd.DataFrame(results)
    output_df.to_csv('master_geospatial_metadata.csv', index=False)
    
    # 🛡️ THE PURGE: Keep ONLY strains mathematically anchored to India
    indian_df = output_df[output_df['Exact_Location'].str.contains('India', na=False, case=False)]
    indian_df.to_csv('pure_indian_metadata.csv', index=False)
    
    print(f"\n[COMPLETE] Interrogation finished.")
    print(f"Total strains processed: {len(output_df)}")
    print(f"Verified Indian clinical strains secured: {len(indian_df)}")

if __name__ == "__main__":
    main()
