#!/usr/bin/env python3

from __future__ import annotations
import re
import time
import zipfile
from pathlib import Path
import pandas as pd
import requests
SCRIPT_FOLDER = Path(__file__).resolve().parent
DATA_FOLDER = SCRIPT_FOLDER / 'data'
INPUT_FILE = DATA_FOLDER / 'Supplementary File 1.xlsx'
SHEET_NAME = 'Table S5'
ENSEMBL_COLUMN = 'Ensembl_ID'
INTACT_ZIP = DATA_FOLDER / 'intact_human.zip'
MAPPING_FILE = DATA_FOLDER / 'ensembl_uniprot_mapping.csv'
OUTPUT_EDGES = DATA_FOLDER / 'intact_ppi_edges.csv'
OUTPUT_UNMAPPED = DATA_FOLDER / 'unmapped_ensembl_ids.csv'
OUTPUT_SCORE_SUMMARY = DATA_FOLDER / 'intact_edge_counts_by_score.csv'
SCAN_MIN_MI_SCORE = 0.1
SCAN_MAX_MI_SCORE = 1.0
FINAL_MIN_MI_SCORE = 0.6
FINAL_MAX_MI_SCORE = 1.0
HUMAN_TAXID = '9606'
INTACT_URL = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip'
UNIPROT_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
ENSEMBL_LOOKUP_URL = 'https://rest.ensembl.org/lookup/symbol/'
REQUEST_TIMEOUT = 90
MAX_RETRIES = 5
CHUNK_SIZE = 100000
session = requests.Session()
session.headers.update({'User-Agent': 'IntAct-local-filter-script/1.0'})

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

def map_one_ensembl_gene(ensembl_id):
    params = {
        'query': f'xref:ensembl-{ensembl_id} AND organism_id:9606',
        'format': 'tsv',
        'fields': 'accession,gene_primary,reviewed',
        'size': 100
    }
    try:
        response = request_with_retry(UNIPROT_SEARCH_URL, params=params)
    except Exception as exc:
        print(f'  Mapping error {ensembl_id}: {exc}')
        return []
    lines = response.text.strip().splitlines()
    if len(lines) <= 1:
        return []
    results = []
    for line in lines[1:]:
        fields = line.split('\t')
        accession = fields[0].strip() if len(fields) > 0 else ''
        label = fields[1].strip() if len(fields) > 1 else ''
        reviewed = fields[2].strip() if len(fields) > 2 else ''
        if not accession:
            continue
        results.append({'Ensembl_ID': ensembl_id, 'UniProt': accession, 'Gene_Label': label, 'Reviewed': reviewed})
    return results

def build_mapping(ensembl_ids):
    existing = pd.DataFrame()
    if MAPPING_FILE.exists():
        print('\nExisting mapping cache found:')
        print(MAPPING_FILE)
        existing = pd.read_csv(MAPPING_FILE, dtype=str).fillna('')
    already_done = set()
    if not existing.empty:
        already_done = set(existing['Ensembl_ID'])
    remaining = [x for x in ensembl_ids if x not in already_done]
    print(f'\nEnsembl IDs already mapped: {len(already_done):,}')
    print(f'Ensembl IDs still to query: {len(remaining):,}')
    new_rows = []
    unmapped = []
    for i, ensembl_id in enumerate(remaining, start=1):
        print(f'[mapping {i}/{len(remaining)}] {ensembl_id}')
        rows = map_one_ensembl_gene(ensembl_id)
        if rows:
            print('    ' + ', '.join((r['UniProt'] for r in rows[:5])))
            new_rows.extend(rows)
        else:
            print('    no UniProt mapping')
            unmapped.append(ensembl_id)
        time.sleep(0.05)
        if len(new_rows) > 0 and i % 100 == 0:
            temp_new = pd.DataFrame(new_rows)
            combined = pd.concat([existing, temp_new], ignore_index=True)
            combined = combined.drop_duplicates(subset=['Ensembl_ID', 'UniProt'])
            combined.to_csv(MAPPING_FILE, index=False)
    if new_rows:
        new_df = pd.DataFrame(new_rows)
        mapping = pd.concat([existing, new_df], ignore_index=True)
    else:
        mapping = existing.copy()
    if not mapping.empty:
        mapping = mapping.drop_duplicates(subset=['Ensembl_ID', 'UniProt']).reset_index(drop=True)
        mapping.to_csv(MAPPING_FILE, index=False)
    if unmapped:
        pd.DataFrame({'Ensembl_ID': unmapped}).to_csv(OUTPUT_UNMAPPED, index=False)
    print(f'\nTotal Ensembl-UniProt mappings: {len(mapping):,}')
    print(f"Unique mapped Ensembl genes: {mapping['Ensembl_ID'].nunique():,}")
    return mapping

def download_intact():
    DATA_FOLDER.mkdir(parents=True, exist_ok=True)
    if INTACT_ZIP.exists():
        size_gb = INTACT_ZIP.stat().st_size / 1024 ** 3
        print(f'\nIntAct archive already exists:\n{INTACT_ZIP}\nSize: {size_gb:.2f} GB')
        return
    print('\nDownloading IntAct human MITAB archive...')
    print(INTACT_URL)
    response = request_with_retry(INTACT_URL, stream=True)
    total = int(response.headers.get('content-length', 0))
    downloaded = 0
    with open(INTACT_ZIP, 'wb') as output:
        for block in response.iter_content(chunk_size=1024 * 1024):
            if not block:
                continue
            output.write(block)
            downloaded += len(block)
            if total:
                pct = downloaded / total * 100
                print(f'\rDownloaded {downloaded / 1024 ** 3:.2f} GB ({pct:.1f}%)', end='')
            else:
                print(f'\rDownloaded {downloaded / 1024 ** 3:.2f} GB', end='')
    print('\nDownload complete.')

def normalize_uniprot(value):
    if pd.isna(value):
        return None
    value = str(value)
    match = re.search('uniprotkb:([A-Z0-9]+(?:-\\d+)?)', value, re.I)
    if not match:
        return None
    accession = match.group(1).upper()
    accession = re.sub('-\\d+$', '', accession)
    return accession

def extract_mi_score(value):
    if pd.isna(value):
        return None
    match = re.search('intact-miscore:([0-9]*\\.?[0-9]+)', str(value), re.I)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None

def extract_gene_label(value):
    if pd.isna(value):
        return ''
    text = str(value)
    patterns = [
        'uniprotkb:([^|()]+)\\(gene name\\)',
        'uniprotkb:([^|()]+)\\(display_short\\)',
        ':([^|()]+)\\(gene name\\)'
    ]
    for pattern in patterns:
        match = re.search(pattern, text, re.I)
        if match:
            return match.group(1).strip()
    return ''

def is_human_taxid(value):
    if pd.isna(value):
        return False
    return bool(re.search('taxid:9606(?:\\D|$)', str(value), re.I))

def is_protein_type(value):
    if pd.isna(value):
        return True
    value = str(value).lower()
    if value in {'', '-', 'nan'}:
        return True
    return 'mi:0326' in value or 'protein' in value

def is_experimental_method(value):
    if pd.isna(value):
        return False
    value = str(value).lower()
    if value in {'', '-', 'nan'}:
        return False
    excluded = ['mi:0363', 'inferred by curator', 'predicted', 'prediction', 'inferred interaction']
    return not any((x in value for x in excluded))

def clean_reference(value):
    if pd.isna(value):
        return ''
    values = []
    for item in str(value).split('|'):
        item = item.strip()
        if item and item != '-' and (item not in values):
            values.append(item)
    return ';'.join(values)

def filter_intact(mapping):
    mapping = mapping.copy()
    mapping['UniProt_norm'] = mapping['UniProt'].astype(str).str.replace('-\\d+$', '', regex=True).str.upper()
    query_uniprots = set(mapping['UniProt_norm'])
    print(f'\nUnique UniProt proteins in query set: {len(query_uniprots):,}')
    uniprot_to_label = mapping.drop_duplicates('UniProt_norm').set_index('UniProt_norm')['Gene_Label'].to_dict()
    uniprot_to_ensembl = mapping.drop_duplicates('UniProt_norm').set_index('UniProt_norm')['Ensembl_ID'].to_dict()
    required_columns = list(range(22))
    output_rows = []
    total_rows = 0
    score_pass = 0
    human_pass = 0
    protein_pass = 0
    experimental_pass = 0
    query_pass = 0
    with zipfile.ZipFile(INTACT_ZIP, 'r') as zf:
        names = zf.namelist()
        txt_files = [n for n in names if n.lower().endswith('.txt')]
        if not txt_files:
            raise RuntimeError('No .txt file found inside IntAct ZIP.')
        print('\nFiles inside IntAct archive:')
        for name in txt_files:
            print(f'  {name}')
        mitab_member = txt_files[0]
        print(f'\nReading:\n{mitab_member}')
        with zf.open(mitab_member, 'r') as raw_file:
            reader = pd.read_csv(
                raw_file,
                sep='\t',
                header=None,
                usecols=required_columns,
                dtype=str,
                chunksize=CHUNK_SIZE,
                low_memory=False,
                on_bad_lines='skip'
            )
            for chunk_number, chunk in enumerate(reader, start=1):
                total_rows += len(chunk)
                human_mask = chunk[9].fillna('').str.contains('taxid:9606(?:\\D|$)', regex=True, case=False) & chunk[10].fillna('').str.contains('taxid:9606(?:\\D|$)', regex=True, case=False)
                chunk = chunk[human_mask].copy()
                human_pass += len(chunk)
                if chunk.empty:
                    continue
                chunk['MI_score'] = chunk[14].str.extract('intact-miscore:([0-9]*\\.?[0-9]+)', expand=False)
                chunk['MI_score'] = pd.to_numeric(chunk['MI_score'], errors='coerce')
                score_mask = chunk['MI_score'].between(SCAN_MIN_MI_SCORE, SCAN_MAX_MI_SCORE, inclusive='both')
                chunk = chunk[score_mask].copy()
                score_pass += len(chunk)
                if chunk.empty:
                    continue
                method = chunk[6].fillna('').str.lower()
                experimental_mask = ~method.str.contains(
                    'mi:0363|inferred by curator|predicted|prediction|inferred interaction',
                    regex=True
                ) & method.ne('') & method.ne('-')
                chunk = chunk[experimental_mask].copy()
                experimental_pass += len(chunk)
                if chunk.empty:
                    continue
                type_a = chunk[20].fillna('').str.lower()
                type_b = chunk[21].fillna('').str.lower()
                protein_a = type_a.eq('') | type_a.eq('-') | type_a.str.contains('mi:0326|protein', regex=True)
                protein_b = type_b.eq('') | type_b.eq('-') | type_b.str.contains('mi:0326|protein', regex=True)
                chunk = chunk[protein_a & protein_b].copy()
                protein_pass += len(chunk)
                if chunk.empty:
                    continue
                combined_a = chunk[0].fillna('') + '|' + chunk[2].fillna('')
                combined_b = chunk[1].fillna('') + '|' + chunk[3].fillna('')
                chunk['UniProt_A'] = combined_a.str.extract('uniprotkb:([A-Z0-9]+(?:-\\d+)?)', expand=False).str.upper().str.replace('-\\d+$', '', regex=True)
                chunk['UniProt_B'] = combined_b.str.extract('uniprotkb:([A-Z0-9]+(?:-\\d+)?)', expand=False).str.upper().str.replace('-\\d+$', '', regex=True)
                input_a = chunk['UniProt_A'].isin(query_uniprots)
                input_b = chunk['UniProt_B'].isin(query_uniprots)
                chunk = chunk[input_a | input_b].copy()
                query_pass += len(chunk)
                if chunk.empty:
                    print(
                        f'\rChunk {chunk_number:,} | rows scanned {total_rows:,} | matches {query_pass:,}',
                        end=''
                    )
                    continue
                for row in chunk.itertuples(index=False, name=None):
                    id_a = row[0]
                    id_b = row[1]
                    alt_a = row[2]
                    alt_b = row[3]
                    alias_a = row[4]
                    alias_b = row[5]
                    detection_method = row[6]
                    first_author = row[7]
                    publication = row[8]
                    interaction_type = row[11]
                    source_database = row[12]
                    interaction_id = row[13]
                    score = row[22]
                    uniprot_a = row[23]
                    uniprot_b = row[24]
                    a_is_input = uniprot_a in query_uniprots
                    b_is_input = uniprot_b in query_uniprots
                    if a_is_input:
                        source_uniprot = uniprot_a
                        target_uniprot = uniprot_b
                        source_alias = alias_a
                        target_alias = alias_b
                    else:
                        source_uniprot = uniprot_b
                        target_uniprot = uniprot_a
                        source_alias = alias_b
                        target_alias = alias_a
                    source_label = uniprot_to_label.get(source_uniprot, '')
                    if not source_label:
                        source_label = extract_gene_label(source_alias)
                    source_ensembl = uniprot_to_ensembl.get(source_uniprot, '')
                    target_label = extract_gene_label(target_alias)
                    output_rows.append({'source_label': source_label, 'source_uniprot': source_uniprot, 'source_ensembl': source_ensembl, 'target_uniprot': target_uniprot, 'target_label': target_label, 'intact_miscore': score, 'detection_method': detection_method, 'interaction_type': interaction_type, 'first_author': first_author, 'publication': clean_reference(publication), 'intact_interaction_id': clean_reference(interaction_id), 'source_database': source_database})
                print(
                    f'\rChunk {chunk_number:,} | rows scanned {total_rows:,} | matches {query_pass:,}',
                    end=''
                )
    print()
    print('\nFiltering summary')
    print(f'Rows scanned:             {total_rows:,}')
    print(f'Human-human:              {human_pass:,}')
    print(f'MI score >= {SCAN_MIN_MI_SCORE:.2f}:        {score_pass:,}')
    print(f'Experimental evidence:    {experimental_pass:,}')
    print(f'Protein-protein:          {protein_pass:,}')
    print(f'Touches input gene list:  {query_pass:,}')
    return pd.DataFrame(output_rows)

def map_partner_labels(evidence):
    if evidence.empty:
        return evidence
    target_labels = evidence['target_label'].fillna('').astype(str).str.strip().replace('', pd.NA).dropna().drop_duplicates().tolist()
    print(f'\nUnique target gene labels to map: {len(target_labels):,}')
    mapping = {}
    for i, label in enumerate(target_labels, start=1):
        print(f'\rMapping target {i:,}/{len(target_labels):,}', end='')
        url = ENSEMBL_LOOKUP_URL + f'homo_sapiens/{label}'
        try:
            response = session.get(
                url,
                headers={'Content-Type': 'application/json', 'Accept': 'application/json'},
                timeout=REQUEST_TIMEOUT
            )
            response.raise_for_status()
            result = response.json()
            ensembl_id = result.get('id', '')
            mapping[label] = ensembl_id
        except Exception:
            mapping[label] = ''
        time.sleep(0.03)
    print()
    evidence['target_ensembl'] = evidence['target_label'].fillna('').astype(str).str.strip().map(mapping).fillna('')
    return evidence

def join_unique(series):
    values = []
    for item in series.dropna():
        for value in str(item).split(';'):
            value = value.strip()
            if value and value != '-' and (value not in values):
                values.append(value)
    return ';'.join(values)

def count_unique(series):
    values = set()
    for item in series.dropna():
        for value in str(item).split(';'):
            value = value.strip()
            if value and value != '-':
                values.add(value)
    return len(values)

def collapse_edges(evidence):
    if evidence.empty:
        return evidence
    evidence = evidence.copy()
    evidence['_node1'] = evidence.apply(lambda r: min(str(r['source_uniprot']), str(r['target_uniprot'])), axis=1)
    evidence['_node2'] = evidence.apply(lambda r: max(str(r['source_uniprot']), str(r['target_uniprot'])), axis=1)
    collapsed = evidence.groupby(
        ['_node1', '_node2'],
        as_index=False
    ).agg(source_label=('source_label', join_unique), source_uniprot=('source_uniprot', join_unique), source_ensembl=('source_ensembl', join_unique), target_ensembl=('target_ensembl', join_unique), target_label=('target_label', join_unique), target_uniprot=('target_uniprot', join_unique), intact_miscore=('intact_miscore', 'max'), evidence_count=('intact_interaction_id', 'size'), detection_method=('detection_method', join_unique), interaction_type=('interaction_type', join_unique), publication=('publication', join_unique), publication_count=('publication', count_unique), intact_interaction_id=('intact_interaction_id', join_unique), first_author=('first_author', join_unique), source_database=('source_database', join_unique))
    collapsed = collapsed.drop(columns=['_node1', '_node2'])
    return collapsed
SCORE_THRESHOLDS = [0.1, 0.4, 0.7, 0.9]

def count_edges_by_score(edges, thresholds=SCORE_THRESHOLDS):
    if edges.empty:
        print('\nNo edges to summarize by MI score.')
        return pd.DataFrame(columns=['min_mi_score', 'edge_count'])
    scores = pd.to_numeric(edges['intact_miscore'], errors='coerce')
    print('\nEdge counts by MI score threshold')
    rows = []
    for threshold in thresholds:
        count = int((scores >= threshold).sum())
        print(f'  MI score >= {threshold:.2f}: {count:,}')
        rows.append({'min_mi_score': threshold, 'edge_count': count})
    return pd.DataFrame(rows)

def main():
    DATA_FOLDER.mkdir(parents=True, exist_ok=True)
    print('========================================')
    print('IntAct human PPI local filtering')
    print('========================================')
    ensembl_ids = read_input_genes()
    mapping = build_mapping(ensembl_ids)
    if mapping.empty:
        raise RuntimeError('No Ensembl IDs could be mapped to UniProt.')
    download_intact()
    evidence = filter_intact(mapping)
    if evidence.empty:
        print('\nNo interactions passed the filters.')
        return
    evidence = map_partner_labels(evidence)
    edges = collapse_edges(evidence)
    edges = edges.sort_values('intact_miscore', ascending=False).reset_index(drop=True)
    score_summary = count_edges_by_score(edges)
    score_summary.to_csv(OUTPUT_SCORE_SUMMARY, index=False)
    edges = edges[edges['intact_miscore'].between(FINAL_MIN_MI_SCORE, FINAL_MAX_MI_SCORE, inclusive='both')].reset_index(drop=True)
    edges.to_csv(OUTPUT_EDGES, index=False)
    print(f'\nEdges kept after final MI score cutoff ({FINAL_MIN_MI_SCORE:.2f}-{FINAL_MAX_MI_SCORE:.2f}): {len(edges):,}')
    print('\n========================================')
    print('FINISHED')
    print('========================================')
    print(f'Unique PPI edges: {len(edges):,}')
    print(f"Input proteins represented: {evidence['source_uniprot'].nunique():,}")
    mapped_targets = edges['target_ensembl'].fillna('').ne('').sum()
    total_targets = len(edges)
    print(f'Edges with target Ensembl ID: {mapped_targets:,} / {total_targets:,}')
    print('\nUnique edge file:')
    print(OUTPUT_EDGES)
    print('\nEdge counts by MI score:')
    print(OUTPUT_SCORE_SUMMARY)
    print('\nIntAct archive:')
    print(INTACT_ZIP)
if __name__ == '__main__':
    main()
