#!/usr/bin/env python3
"""
NCBI Genome Extractor - Integrated version for MutationScan
"""

import argparse
import csv
import json
import logging
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import List, Optional, Dict, Any
from urllib.parse import parse_qs, urlparse

import requests
from tqdm import tqdm

# Configuration
NCBI_EMAIL = "vihaankulkarni29@gmail.com"
NCBI_API_KEY = "ef7622c2e716fa317fe04d24c42904211107"
DEFAULT_MAX_RESULTS = 100
DEFAULT_RETRIES = 3
DEFAULT_DELAY = 0.5
DEFAULT_OUTPUT_DIR = "./output"
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_SEARCH_URL = NCBI_BASE_URL + "esearch.fcgi"
NCBI_FETCH_URL = NCBI_BASE_URL + "efetch.fcgi"


class NCBIExtractor:
    """Main class for NCBI genome extraction - integrated for MutationScan"""

    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY):
        self.email = email
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': f'MutationScan-NCBIExtractor/1.0 (mailto:{self.email})'
        })

    def search_genomes(self, query: str, max_results: int = DEFAULT_MAX_RESULTS,
                      select_best: bool = True) -> List[str]:
        """Search NCBI for genome IDs matching the query"""
        if select_best:
            initial_results = max_results * 5
            initial_results = min(initial_results, 500)
        else:
            initial_results = max_results

        params = {
            'db': 'nuccore',
            'term': query,
            'retmax': initial_results,
            'retmode': 'xml',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}
        response = self._make_request(NCBI_SEARCH_URL, params)

        if not response:
            return []

        root = ET.fromstring(response.text)
        ids = []
        for id_element in root.findall('.//Id'):
            ids.append(id_element.text)

        if select_best and len(ids) > max_results:
            ids = self._select_best_genomes(ids, max_results)

        return ids

    def download_fasta(self, genome_id: str, output_dir: str = DEFAULT_OUTPUT_DIR,
                      retries: int = DEFAULT_RETRIES) -> bool:
        """Download FASTA sequence for a single genome ID"""
        os.makedirs(output_dir, exist_ok=True)

        params = {
            'db': 'nuccore',
            'id': genome_id,
            'rettype': 'fasta',
            'retmode': 'text',
            'email': self.email,
            'api_key': self.api_key if self.api_key else None
        }

        params = {k: v for k, v in params.items() if v is not None}

        for attempt in range(retries + 1):
            try:
                response = self._make_request(NCBI_FETCH_URL, params, stream=True)
                if not response:
                    continue

                filename = f"{genome_id}.fasta"
                filepath = os.path.join(output_dir, filename)

                total_size = int(response.headers.get('content-length', 0))
                with open(filepath, 'wb') as f, tqdm(
                    desc=f"Downloading {genome_id}",
                    total=total_size,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024,
                ) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

                logging.info(f"Successfully downloaded {genome_id}")
                return True

            except Exception as e:
                logging.warning(f"Attempt {attempt + 1} failed for {genome_id}: {e}")
                if attempt < retries:
                    time.sleep(DEFAULT_DELAY * (2 ** attempt))

        logging.error(f"Failed to download {genome_id} after {retries + 1} attempts")
        return False

    def download_multiple_genomes(self, genome_ids: List[str], output_dir: str = DEFAULT_OUTPUT_DIR,
                                max_concurrent: int = 3) -> List[str]:
        """Download multiple genomes with progress tracking"""
        successful = []
        failed = []

        with tqdm(total=len(genome_ids), desc="Overall Progress") as pbar:
            for genome_id in genome_ids:
                if self.download_fasta(genome_id, output_dir):
                    successful.append(genome_id)
                else:
                    failed.append(genome_id)
                pbar.update(1)

        logging.info(f"Downloaded {len(successful)} genomes successfully")
        if failed:
            logging.warning(f"Failed to download {len(failed)} genomes: {failed}")

        return successful

    def _make_request(self, url: str, params: dict, stream: bool = False) -> Optional[requests.Response]:
        """Make HTTP request with error handling"""
        try:
            response = self.session.get(url, params=params, stream=stream, timeout=30)
            response.raise_for_status()
            return response
        except requests.RequestException as e:
            logging.error(f"Request failed: {e}")
            return None

    def _select_best_genomes(self, genome_ids: List[str], max_results: int) -> List[str]:
        """Select the best genomes based on metadata quality"""
        logging.info(f"Evaluating {len(genome_ids)} genomes for metadata quality...")
        # For now, return first N - can be enhanced later
        return genome_ids[:max_results]


# Main function for backward compatibility
def extract_genomes(query: str, output_dir: str, max_results: int = 10) -> List[str]:
    """Simple function interface for genome extraction"""
    extractor = NCBIExtractor()
    genome_ids = extractor.search_genomes(query, max_results)
    if genome_ids:
        return extractor.download_multiple_genomes(genome_ids, output_dir)
    return []


if __name__ == '__main__':
    # Simple CLI interface
    import argparse
    parser = argparse.ArgumentParser(description="Download genomes from NCBI")
    parser.add_argument('--query', required=True, help='NCBI search query')
    parser.add_argument('--output-dir', default='./genomes', help='Output directory')
    parser.add_argument('--max-results', type=int, default=10, help='Maximum genomes to download')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    extract_genomes(args.query, args.output_dir, args.max_results)