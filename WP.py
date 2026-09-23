#!/usr/bin/env python3

from __future__ import annotations
import re
import time
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import pandas as pd
import requests
SCRIPT_FOLDER = Path(__file__).resolve().parent
DATA_FOLDER = SCRIPT_FOLDER / 'v3data'
INPUT_FILE = DATA_FOLDER / 'OMIM_gene_mapped_comparison.xlsx'
SHEET_NAME = 'Summary'
ENSEMBL_COLUMN = 'Ensembl_ID'
WIKIPATHWAYS_ZIP = DATA_FOLDER / 'wikipathways_human_gpml.zip'
OUTPUT_EDGES = DATA_FOLDER / 'wikipathways_ppi_edges.csv'
OUTPUT_NOT_IN_WIKIPATHWAYS = DATA_FOLDER / 'genes_not_in_wikipathways.csv'
HUMAN_ORGANISM = 'Homo sapiens'
WIKIPATHWAYS_GPML_INDEX_URL = 'https://data.wikipathways.org/current/gpml/'
ALLOWED_DATANODE_TYPES = {'geneproduct', 'protein'}
PROTEIN_LEVEL_INTERACTION_TYPES = {
    'mim-binding',
    'mim-complex',
    'mim-catalysis',
    'mim-stimulation',
    'mim-inhibition',
    'mim-necessary-stimulation',
    'mim-modification'
}
REQUEST_TIMEOUT = 90
MAX_RETRIES = 5
PROGRESS_EVERY = 200
session = requests.Session()
session.headers.update({'User-Agent': 'WikiPathways-local-filter-script/1.0'})

def request_with_retry(url, params=None, stream=False):
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            response = session.get(url, params=params, timeout=REQUEST_TIMEOUT, stream=stream)
            if response.status_code == 429:
                wait = 2 ** (attempt + 1)
                print(f'Rate limited. Waiting {wait} seconds...')
                time.sleep(wait)
                continue
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
            if attempt == MAX_RETRIES - 1:
                break
            wait = 2 ** attempt
            print(f'Request failed: {exc}\nRetrying in {wait}s...')
            time.sleep(wait)
    raise RuntimeError(f'Request failed after retries: {last_error}')

def download_file(url, dest_path):
    if dest_path.exists():
        size_mb = dest_path.stat().st_size / 1024 ** 2
        print(f'\nAlready downloaded:\n{dest_path}\nSize: {size_mb:.1f} MB')
        return
    print(f'\nDownloading:\n{url}')
    response = request_with_retry(url, stream=True)
    total = int(response.headers.get('content-length', 0))
    downloaded = 0
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, 'wb') as output:
        for block in response.iter_content(chunk_size=1024 * 1024):
            if not block:
                continue
            output.write(block)
            downloaded += len(block)
            if total:
                pct = downloaded / total * 100
                print(f'\rDownloaded {downloaded / 1024 ** 2:.1f} MB ({pct:.1f}%)', end='')
            else:
                print(f'\rDownloaded {downloaded / 1024 ** 2:.1f} MB', end='')
    print('\nDownload complete.')

def read_input_genes():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f'Input file not found:\n{INPUT_FILE}')
    print(f'Reading:\n{INPUT_FILE}')
    df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME)
    df.columns = [str(c).strip() for c in df.columns]
    if ENSEMBL_COLUMN not in df.columns:
        raise ValueError(f"Column '{ENSEMBL_COLUMN}' not found.\n\nColumns found:\n{list(df.columns)}")
    ids = df[ENSEMBL_COLUMN].dropna().astype(str).str.strip().str.replace('\\.\\d+$', '', regex=True)
    ids = ids[ids.str.match('^ENSG\\d+$', na=False)]
    ids = ids.drop_duplicates().tolist()
    print(f'Unique valid Ensembl gene IDs: {len(ids):,}')
    return ids

def find_latest_gpml_zip_url():
    print(f'\nLooking up current GPML archive at:\n{WIKIPATHWAYS_GPML_INDEX_URL}')
    response = request_with_retry(WIKIPATHWAYS_GPML_INDEX_URL)
    matches = re.findall(
        'href="(wikipathways-\\d{8}-gpml-Homo_sapiens\\.zip)"',
        response.text,
        re.IGNORECASE
    )
    if not matches:
        raise RuntimeError(f'Could not find a Homo_sapiens GPML archive link on {WIKIPATHWAYS_GPML_INDEX_URL}. WikiPathways may have changed its directory layout - check the URL manually.')
    filename = sorted(matches)[-1]
    url = WIKIPATHWAYS_GPML_INDEX_URL + filename
    print(f'Found archive: {filename}')
    return (url, filename)

def get_namespace(root):
    if root.tag.startswith('{'):
        return root.tag[:root.tag.index('}') + 1]
    return ''

def extract_pathway_id(root, ns, filename):
    xref = root.find(f'{ns}Xref')
    if xref is not None:
        database = (xref.get('Database') or '').strip().lower()
        if 'wikipathways' in database:
            wp_id = xref.get('ID', '').strip()
            if wp_id:
                return wp_id
    match = re.search('(WP\\d+)', filename)
    return match.group(1) if match else filename

def parse_gpml(xml_bytes, filename):
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return []
    ns = get_namespace(root)
    pathway_name = root.get('Name', '').strip()
    pathway_id = extract_pathway_id(root, ns, filename)
    graphid_to_node = {}
    for datanode in root.iter(f'{ns}DataNode'):
        graph_id = datanode.get('GraphId')
        if not graph_id:
            continue
        node_type = (datanode.get('Type') or '').strip().lower()
        if node_type not in ALLOWED_DATANODE_TYPES:
            continue
        label = (datanode.get('TextLabel') or '').strip()
        ensembl_id = ''
        xref = datanode.find(f'{ns}Xref')
        if xref is not None:
            database = (xref.get('Database') or '').strip().lower()
            xref_id = (xref.get('ID') or '').strip()
            if database == 'ensembl' and xref_id:
                ensembl_id = xref_id
        if not ensembl_id:
            continue
        graphid_to_node[graph_id] = {'ensembl': ensembl_id, 'label': label}
    edges = []
    for interaction in root.iter(f'{ns}Interaction'):
        graphics = interaction.find(f'{ns}Graphics')
        if graphics is None:
            continue
        points = graphics.findall(f'{ns}Point')
        refs = [(p.get('GraphRef'), p.get('ArrowHead')) for p in points if p.get('GraphRef')]
        if len(refs) < 2:
            continue
        source_ref = refs[0][0]
        target_ref, arrowhead = (refs[-1][0], refs[-1][1])
        if not source_ref or not target_ref or source_ref == target_ref:
            continue
        source_node = graphid_to_node.get(source_ref)
        target_node = graphid_to_node.get(target_ref)
        if source_node is None or target_node is None:
            continue
        interaction_type = (arrowhead or '').strip().lower() or 'unspecified'
        edges.append({'node_a_ensembl': source_node['ensembl'], 'node_a_label': source_node['label'], 'node_b_ensembl': target_node['ensembl'], 'node_b_label': target_node['label'], 'interaction_type': interaction_type, 'pathway_id': pathway_id, 'pathway_name': pathway_name})
    return edges

def filter_wikipathways(query_ensembl_ids):
    output_rows = []
    seen_ensembl_ids = set()
    with zipfile.ZipFile(WIKIPATHWAYS_ZIP, 'r') as zf:
        gpml_names = [n for n in zf.namelist() if n.lower().endswith('.gpml')]
        print(f'\nPathway files in archive: {len(gpml_names):,}')
        total_edges_seen = 0
        total_edges_kept = 0
        for i, name in enumerate(gpml_names, start=1):
            xml_bytes = zf.read(name)
            edges = parse_gpml(xml_bytes, name)
            total_edges_seen += len(edges)
            for edge in edges:
                a_ensembl = edge['node_a_ensembl']
                b_ensembl = edge['node_b_ensembl']
                seen_ensembl_ids.add(a_ensembl)
                seen_ensembl_ids.add(b_ensembl)
                a_is_input = a_ensembl in query_ensembl_ids
                b_is_input = b_ensembl in query_ensembl_ids
                if not (a_is_input or b_is_input):
                    continue
                if PROTEIN_LEVEL_INTERACTION_TYPES is not None and edge['interaction_type'] not in PROTEIN_LEVEL_INTERACTION_TYPES:
                    continue
                if a_is_input:
                    source_ensembl, source_label = (a_ensembl, edge['node_a_label'])
                    target_ensembl, target_label = (b_ensembl, edge['node_b_label'])
                else:
                    source_ensembl, source_label = (b_ensembl, edge['node_b_label'])
                    target_ensembl, target_label = (a_ensembl, edge['node_a_label'])
                output_rows.append({'source_label': source_label, 'source_ensembl': source_ensembl, 'target_label': target_label, 'target_ensembl': target_ensembl, 'interaction_type': edge['interaction_type'], 'pathway_id': edge['pathway_id'], 'pathway_name': edge['pathway_name'], 'source_database': 'WikiPathways'})
                total_edges_kept += 1
            if i % PROGRESS_EVERY == 0 or i == len(gpml_names):
                print(
                    f'\rPathways processed: {i:,}/{len(gpml_names):,} | edges seen {total_edges_seen:,} | edges kept {total_edges_kept:,}',
                    end=''
                )
    print()
    print('\nFiltering summary')
    print(f'Pathway files processed:   {len(gpml_names):,}')
    print(f'Total interactions seen:   {total_edges_seen:,}')
    print(f'Interactions kept:         {total_edges_kept:,}')
    not_in_wikipathways = [x for x in query_ensembl_ids if x not in seen_ensembl_ids]
    if not_in_wikipathways:
        pd.DataFrame({'Ensembl_ID': not_in_wikipathways}).to_csv(OUTPUT_NOT_IN_WIKIPATHWAYS, index=False)
    print(f'Input genes never seen as a WikiPathways DataNode: {len(not_in_wikipathways):,}')
    result = pd.DataFrame(output_rows)
    if not result.empty:
        result = result[result['source_ensembl'].astype(str).str.strip().ne('') & result['target_ensembl'].astype(str).str.strip().ne('')].reset_index(drop=True)
    return result

def join_unique(series):
    values = []
    for item in series.dropna():
        for value in str(item).split(';'):
            value = value.strip()
            if value and value not in values:
                values.append(value)
    return ';'.join(values)

def count_unique(series):
    values = set()
    for item in series.dropna():
        for value in str(item).split(';'):
            value = value.strip()
            if value:
                values.add(value)
    return len(values)

def collapse_edges(evidence):
    if evidence.empty:
        return evidence
    evidence = evidence.copy()
    evidence['_node1'] = evidence[['source_ensembl', 'target_ensembl']].min(axis=1)
    evidence['_node2'] = evidence[['source_ensembl', 'target_ensembl']].max(axis=1)
    collapsed = evidence.groupby(
        ['_node1', '_node2'],
        as_index=False
    ).agg(source_label=('source_label', join_unique), source_ensembl=('source_ensembl', join_unique), target_label=('target_label', join_unique), target_ensembl=('target_ensembl', join_unique), interaction_type=('interaction_type', join_unique), pathway_id=('pathway_id', join_unique), pathway_count=('pathway_id', count_unique), pathway_name=('pathway_name', join_unique), evidence_count=('pathway_id', 'size'), source_database=('source_database', join_unique))
    collapsed = collapsed.drop(columns=['_node1', '_node2'], errors='ignore')
    return collapsed

def main():
    DATA_FOLDER.mkdir(parents=True, exist_ok=True)
    print('=' * 40)
    print('WikiPathways human PPI local filtering')
    print('=' * 40)
    ensembl_ids = read_input_genes()
    query_ensembl_ids = set(ensembl_ids)
    zip_url, zip_filename = find_latest_gpml_zip_url()
    download_file(zip_url, WIKIPATHWAYS_ZIP)
    evidence = filter_wikipathways(query_ensembl_ids)
    if evidence.empty:
        print('\nNo interactions passed the filters.')
        return
    edges = collapse_edges(evidence)
    edges = edges.sort_values('pathway_count', ascending=False).reset_index(drop=True)
    edges.to_csv(OUTPUT_EDGES, index=False)
    print('\n' + '=' * 40)
    print('FINISHED')
    print('=' * 40)
    print(f'Unique PPI edges: {len(edges):,}')
    print('\nUnique edge file:')
    print(OUTPUT_EDGES)
if __name__ == '__main__':
    main()
