#!/usr/bin/env python3

from __future__ import annotations

import gzip
import re
import time
from pathlib import Path
import pandas as pd
import requests


# ============================================================
# Input and Output
# ============================================================

SCRIPT_FOLDER = Path(__file__).resolve().parent
DATA_FOLDER = SCRIPT_FOLDER / "v3data"

# Input Supplementary file for the gene input file
INPUT_FILE = DATA_FOLDER / "Supplementary File 1.xlsx"
SHEET_NAME = "Table S5" 
ENSEMBL_COLUMN = "Ensembl_ID"

# Downloaded STRING files 
STRING_ALIASES_GZ = DATA_FOLDER / "9606.protein.aliases.txt.gz"
STRING_INFO_GZ = DATA_FOLDER / "9606.protein.info.txt.gz"
STRING_LINKS_FULL_GZ = DATA_FOLDER / "9606.protein.links.full.txt.gz"

# Mapping cache: Ensembl gene -> STRING protein id
MAPPING_FILE = DATA_FOLDER / "ensembl_string_mapping.csv"

# Output: unique PPI edges only
OUTPUT_EDGES = DATA_FOLDER / "string_ppi_edges.csv"

# Unmapped input Ensembl IDs
OUTPUT_UNMAPPED = DATA_FOLDER / "unmapped_ensembl_ids_string.csv"

# Edge counts broken down by combined_score threshold
OUTPUT_SCORE_SUMMARY = DATA_FOLDER / "string_edge_counts_by_score.csv"


# ============================================================
# FILTER SETTINGS
# ============================================================

# STRING combined_score is 0-1000.

SCAN_MIN_SCORE = 100
SCAN_MAX_SCORE = 1000

FINAL_MIN_SCORE = 990 ## 990 for high confidence
FINAL_MAX_SCORE = 1000
SCORE_THRESHOLDS = [100, 400, 600, 700, 900, 950, 990] ## To check the number of edges at different confidence levels

REQUIRE_EXPERIMENTAL_EVIDENCE = True## Only experimental evidence

HUMAN_TAXID = "9606" ## Only human PPIs
STRING_VERSION = "12.0"

# STRING bulk download URLs.
STRING_ALIASES_URL = (
    f"https://stringdb-downloads.org/download/"
    f"protein.aliases.v{STRING_VERSION}/"
    f"{HUMAN_TAXID}.protein.aliases.v{STRING_VERSION}.txt.gz"
)

STRING_INFO_URL = (
    f"https://stringdb-downloads.org/download/"
    f"protein.info.v{STRING_VERSION}/"
    f"{HUMAN_TAXID}.protein.info.v{STRING_VERSION}.txt.gz"
)

STRING_LINKS_FULL_URL = (
    f"https://stringdb-downloads.org/download/"
    f"protein.links.full.v{STRING_VERSION}/"
    f"{HUMAN_TAXID}.protein.links.full.v{STRING_VERSION}.txt.gz"
)

REQUEST_TIMEOUT = 90
MAX_RETRIES = 5
CHUNK_SIZE = 200_000
# ============================================================
# HTTP SESSION
# ============================================================

session = requests.Session()
session.headers.update({"User-Agent": "STRING-local-filter-script/1.0"})


def normalize_string_id(value):
    value = str(value).strip()
    if not value:
        return value
    if not value.startswith(f"{HUMAN_TAXID}."):
        value = f"{HUMAN_TAXID}.{value}"
    return value.upper()


def request_with_retry(url, params=None, stream=False):
    last_error = None
    for attempt in range(MAX_RETRIES):
        try:
            response = session.get(
                url, params=params, timeout=REQUEST_TIMEOUT, stream=stream
            )
            if response.status_code == 429:
                wait = 2 ** (attempt + 1)
                print(f"Rate limited. Waiting {wait} seconds...")
                time.sleep(wait)
                continue
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_error = exc
            if attempt == MAX_RETRIES - 1:
                break
            wait = 2**attempt
            print(f"Request failed: {exc}\nRetrying in {wait}s...")
            time.sleep(wait)
    raise RuntimeError(f"Request failed after retries: {last_error}")


def download_file(url, dest_path):
    if dest_path.exists():
        size_mb = dest_path.stat().st_size / 1024**2
        print(f"\nAlready downloaded:\n{dest_path}\nSize: {size_mb:.1f} MB")
        return
    print(f"\nDownloading:\n{url}")
    response = request_with_retry(url, stream=True)
    total = int(response.headers.get("content-length", 0))
    downloaded = 0
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, "wb") as output:
        for block in response.iter_content(chunk_size=1024 * 1024):
            if not block:
                continue
            output.write(block)
            downloaded += len(block)

            if total:
                pct = downloaded / total * 100
                print(
                    f"\rDownloaded {downloaded / 1024**2:.1f} MB "
                    f"({pct:.1f}%)",
                    end="",
                )
            else:
                print(f"\rDownloaded {downloaded / 1024**2:.1f} MB", end="")

    print("\nDownload complete.")

# ============================================================
# READ INPUT ENSEMBL IDS
# ============================================================

def read_input_genes():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found:\n{INPUT_FILE}")
    print(f"Reading:\n{INPUT_FILE}")
    df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME)
    df.columns = [str(c).strip() for c in df.columns]
    if ENSEMBL_COLUMN not in df.columns:
        raise ValueError(
            f"Column '{ENSEMBL_COLUMN}' not found.\n\n"
            f"Columns found:\n{list(df.columns)}"
        )
    ids = (
        df[ENSEMBL_COLUMN]
        .dropna()
        .astype(str)
        .str.strip()
        .str.replace(r"\.\d+$", "", regex=True)
    )
    ids = ids[ids.str.match(r"^ENSG\d+$", na=False)]
    ids = ids.drop_duplicates().tolist()

    print(f"Unique valid Ensembl gene IDs: {len(ids):,}")

    return ids

# ============================================================
# PROTEIN INFO (STRING id -> gene symbol)
# ============================================================

def load_protein_info():
    print(f"\nLoading STRING protein info:\n{STRING_INFO_GZ}")

    info = pd.read_csv(
        STRING_INFO_GZ,
        sep="\t",
        header=None,
        names=["string_id", "preferred_name", "protein_size", "annotation"],
        skiprows=1,
        dtype=str,
        usecols=[0, 1],
    )
    info["string_id"] = info["string_id"].map(normalize_string_id)

    return dict(zip(info["string_id"], info["preferred_name"]))


# ============================================================
# ENSEMBL <-> STRING MAPPING (via aliases file)
# ============================================================

def build_mapping(ensembl_ids):
    if MAPPING_FILE.exists():
        print(f"\nExisting mapping cache found:\n{MAPPING_FILE}")
        mapping = pd.read_csv(MAPPING_FILE, dtype=str).fillna("")
        mapping["STRING_ID"] = mapping["STRING_ID"].map(normalize_string_id)
    else:
        print(f"\nParsing STRING aliases file:\n{STRING_ALIASES_GZ}")
        ensembl_set = set(ensembl_ids)
        rows = []
        reader = pd.read_csv(
            STRING_ALIASES_GZ,
            sep="\t",
            header=None,
            names=["string_id", "alias", "source"],
            skiprows=1,
            dtype=str,
            chunksize=CHUNK_SIZE,
        )
        for i, chunk in enumerate(reader, start=1):
            hits = chunk[chunk["alias"].isin(ensembl_set)]
            if not hits.empty:
                rows.append(
                    hits[["alias", "string_id"]].rename(
                        columns={"alias": "Ensembl_ID", "string_id": "STRING_ID"}
                    )
                )
            print(f"\rAliases file chunks scanned: {i:,}", end="")
        print()
        if rows:
            mapping = pd.concat(rows, ignore_index=True).drop_duplicates()
        else:
            mapping = pd.DataFrame(columns=["Ensembl_ID", "STRING_ID"])
        mapping["STRING_ID"] = mapping["STRING_ID"].map(normalize_string_id)
        mapping = mapping.drop_duplicates()
        mapping.to_csv(MAPPING_FILE, index=False)
    mapped_ids = set(mapping["Ensembl_ID"])
    unmapped = [x for x in ensembl_ids if x not in mapped_ids]
    if unmapped:
        pd.DataFrame({"Ensembl_ID": unmapped}).to_csv(
            OUTPUT_UNMAPPED, index=False
        )
    print(f"\nTotal Ensembl-STRING mappings: {len(mapping):,}")
    print(f"Unique mapped Ensembl genes: {mapping['Ensembl_ID'].nunique():,}")
    print(f"Unmapped Ensembl genes: {len(unmapped):,}")
    return mapping


def build_reverse_ensembl_map():
    print(f"\nBuilding STRING id -> Ensembl gene lookup from:\n{STRING_ALIASES_GZ}")
    reverse_map = {}
    reader = pd.read_csv(
        STRING_ALIASES_GZ,
        sep="\t",
        header=None,
        names=["string_id", "alias", "source"],
        skiprows=1,
        dtype=str,
        chunksize=CHUNK_SIZE,
    )
    ensg_pattern = re.compile(r"^ENSG\d+$")
    for i, chunk in enumerate(reader, start=1):
        hits = chunk[chunk["alias"].str.match(ensg_pattern, na=False)]
        for string_id, ensg in zip(hits["string_id"], hits["alias"]):
            # Keep the first ENSG seen per STRING id
            string_id = normalize_string_id(string_id)
            reverse_map.setdefault(string_id, ensg)
        print(f"\rAliases file chunks scanned: {i:,}", end="")
    print()
    return reverse_map


# ============================================================
# READ AND FILTER STRING LINKS (FULL, WITH EVIDENCE CHANNELS)
# ============================================================

LINKS_COLUMNS = [
    "protein1",
    "protein2",
    "neighborhood",
    "neighborhood_transferred",
    "fusion",
    "cooccurence",
    "homology",
    "coexpression",
    "coexpression_transferred",
    "experimental",
    "experimental_transferred",
    "database",
    "database_transferred",
    "textmining",
    "textmining_transferred",
    "combined_score",
]

LINKS_SCORE_COLUMNS = [c for c in LINKS_COLUMNS if c not in ("protein1", "protein2")]

def filter_string_links(query_string_ids, string_to_label, string_to_ensembl):
    print(f"\nUnique STRING proteins in query set: {len(query_string_ids):,}")
    total_rows = 0
    score_pass = 0
    experimental_pass = 0
    query_pass = 0
    output_rows = []
    dtype_map = {"protein1": str, "protein2": str}
    dtype_map.update({c: "int32" for c in LINKS_SCORE_COLUMNS})
    reader = pd.read_csv(
        STRING_LINKS_FULL_GZ,
        sep=" ",
        header=None,
        names=LINKS_COLUMNS,
        skiprows=1,
        dtype=dtype_map,
        chunksize=CHUNK_SIZE,
    )
    for chunk_number, chunk in enumerate(reader, start=1):
        total_rows += len(chunk)
        chunk["protein1"] = chunk["protein1"].map(normalize_string_id)
        chunk["protein2"] = chunk["protein2"].map(normalize_string_id)
        score_mask = chunk["combined_score"].between(
            SCAN_MIN_SCORE, SCAN_MAX_SCORE, inclusive="both"
        )
        chunk = chunk[score_mask]
        score_pass += len(chunk)
        if chunk.empty:
            continue
        if REQUIRE_EXPERIMENTAL_EVIDENCE:
            chunk = chunk[
                (chunk["experimental"] > 0)
                | (chunk["experimental_transferred"] > 0)
            ]
        experimental_pass += len(chunk)
        if chunk.empty:
            continue
        # ------------------------------------------------
        # At least one participant must be in the input set
        # ------------------------------------------------

        input_1 = chunk["protein1"].isin(query_string_ids)
        input_2 = chunk["protein2"].isin(query_string_ids)
        chunk = chunk[input_1 | input_2]
        query_pass += len(chunk)

        if chunk.empty:
            print(
                f"\rChunk {chunk_number:,} | rows scanned {total_rows:,} "
                f"| matches {query_pass:,}",
                end="",
            )
            continue
        for row in chunk.itertuples(index=False):
            p1 = row.protein1
            p2 = row.protein2

            p1_is_input = p1 in query_string_ids

            if p1_is_input:
                source_id, target_id = p1, p2
            else:
                source_id, target_id = p2, p1
            output_rows.append(
                {
                    "source_label": string_to_label.get(source_id, ""),
                    "source_string_id": source_id,
                    "source_ensembl": string_to_ensembl.get(source_id, ""),
                    "target_string_id": target_id,
                    "target_label": string_to_label.get(target_id, ""),
                    "target_ensembl": string_to_ensembl.get(target_id, ""),
                    "combined_score": row.combined_score,
                    "neighborhood": row.neighborhood,
                    "neighborhood_transferred": row.neighborhood_transferred,
                    "fusion": row.fusion,
                    "cooccurence": row.cooccurence,
                    "homology": row.homology,
                    "coexpression": row.coexpression,
                    "coexpression_transferred": row.coexpression_transferred,
                    "experimental": row.experimental,
                    "experimental_transferred": row.experimental_transferred,
                    "database": row.database,
                    "database_transferred": row.database_transferred,
                    "textmining": row.textmining,
                    "textmining_transferred": row.textmining_transferred,
                    "source_database": "STRING",
                }
            )
        print(
            f"\rChunk {chunk_number:,} | rows scanned {total_rows:,} "
            f"| matches {query_pass:,}",
            end="",
        )

    print()
    print("\nFiltering summary")
    print(f"Rows scanned:             {total_rows:,}")
    print(f"Combined score >= {SCAN_MIN_SCORE}:    {score_pass:,}")
    print(f"Experimental evidence:    {experimental_pass:,}")
    print(f"Touches input gene list:  {query_pass:,}")

    return pd.DataFrame(output_rows)


# ============================================================
# COLLAPSE INTO UNIQUE EDGES
# ============================================================

def collapse_edges(evidence):
    if evidence.empty:
        return evidence
    evidence = evidence.copy()
    evidence["_node1"] = evidence[["source_string_id", "target_string_id"]].min(axis=1)
    evidence["_node2"] = evidence[["source_string_id", "target_string_id"]].max(axis=1)

    collapsed = (
        evidence.sort_values("combined_score", ascending=False)
        .groupby(["_node1", "_node2"], as_index=False)
        .first()
    )
    collapsed = collapsed.drop(columns=["_node1", "_node2"])
    return collapsed


# ============================================================
# COUNT EDGES BY COMBINED SCORE THRESHOLD
# ============================================================

def count_edges_by_score(edges, thresholds=SCORE_THRESHOLDS):
    if edges.empty:
        print("\nNo edges to summarize by combined score.")
        return pd.DataFrame(columns=["min_combined_score", "edge_count"])
    scores = pd.to_numeric(edges["combined_score"], errors="coerce")
    print("\nEdge counts by combined score threshold")
    rows = []
    for threshold in thresholds:
        count = int((scores >= threshold).sum())
        print(f"  Combined score >= {threshold}: {count:,}")
        rows.append({"min_combined_score": threshold, "edge_count": count})
    return pd.DataFrame(rows)


# ============================================================
# MAIN
# ============================================================

def main():
    DATA_FOLDER.mkdir(parents=True, exist_ok=True)

    print("=" * 40)
    print("STRING human PPI local filtering")
    print("=" * 40)
    # 1. Read input ENSG IDs
    ensembl_ids = read_input_genes()

    # 2. Download STRING reference files
    download_file(STRING_ALIASES_URL, STRING_ALIASES_GZ)
    download_file(STRING_INFO_URL, STRING_INFO_GZ)
    download_file(STRING_LINKS_FULL_URL, STRING_LINKS_FULL_GZ)

    # 3. Ensembl -> STRING mapping (for the input gene list)
    mapping = build_mapping(ensembl_ids)
    if mapping.empty:
        raise RuntimeError("No Ensembl IDs could be mapped to STRING proteins.")
    query_string_ids = set(mapping["STRING_ID"])

    # 4. STRING id -> gene symbol (all proteins, for labeling partners)
    string_to_label = load_protein_info()

    # 5. STRING id -> Ensembl gene (all proteins, for labeling partners)
    string_to_ensembl = build_reverse_ensembl_map()

    # 6. Local filtering of the full links file
    evidence = filter_string_links(query_string_ids, string_to_label, string_to_ensembl)
    if evidence.empty:
        print("\nNo interactions passed the filters.")
        return

    # 7. Collapse into unique edges
    edges = collapse_edges(evidence)

    # 8. Sort (full scan-range edges, SCAN_MIN_SCORE-SCAN_MAX_SCORE)
    edges = edges.sort_values("combined_score", ascending=False).reset_index(drop=True)

    # 9. Count edges by combined score threshold.
    # Done on the FULL scan-range edge table so every threshold
    # below shows a real, different count instead of being
    # pre-cut at the final cutoff.
    score_summary = count_edges_by_score(edges)
    score_summary.to_csv(OUTPUT_SCORE_SUMMARY, index=False)

    # 10. Apply the FINAL score cutoff and save only that range
    # as the edge output.
    edges = edges[
        edges["combined_score"].between(
            FINAL_MIN_SCORE, FINAL_MAX_SCORE, inclusive="both"
        )
    ].reset_index(drop=True)
    edges.to_csv(OUTPUT_EDGES, index=False)
    print(
        f"\nEdges kept after final combined score cutoff "
        f"({FINAL_MIN_SCORE}-{FINAL_MAX_SCORE}): {len(edges):,}"
    )
    print("\n" + "=" * 40)
    print("FINISHED")
    print("=" * 40)
    print(f"Unique PPI edges: {len(edges):,}")
    print(f"Input proteins represented: {evidence['source_string_id'].nunique():,}")
    mapped_targets = edges["target_ensembl"].fillna("").ne("").sum()
    print(f"Edges with target Ensembl ID: {mapped_targets:,} / {len(edges):,}")
    print("\nUnique edge file:")
    print(OUTPUT_EDGES)
    print("\nEdge counts by combined score:")
    print(OUTPUT_SCORE_SUMMARY)


if __name__ == "__main__":
    main()
