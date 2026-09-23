# # Genetic Epilepsy Interaction Networks
Scripts for retrieving human gene-associated interactions from STRING, IntAct, and WikiPathways for genetic epilepsy network analysis. The Python scripts take Ensembl gene IDs as input and export filtered, undirected interaction tables for downstream analysis.


## Repository contents

| File | Purpose |
| --- | --- |
| `STRING.py` | Retrieve STRING associations, apply confidence and evidence filters, and export unique protein pairs. |
| `IntAct.py` | Retrieve human IntAct interactions, apply MI-score and detection-method filters, and export unique protein pairs. 
| `WP.py` | Extract selected protein-level interactions from human WikiPathways GPML files. |
| `exclude chemicals.txt` | Reference list of chemicals for exclusion. |
| `__.ipynb` | Notebook scrip for data analysis and visualization. |


## Requirements

The Python scripts require Python 3 and the following packages:

```bash
python -m pip install pandas requests openpyxl
```

## Input data

Place the input workbook in a `v3data` folder beside the scripts. Each script reads one Ensembl gene ID per row from the `Ensembl_ID` column, removes version suffixes, and drops duplicate IDs.

The current default inputs differ between scripts:

| Script | Workbook in `v3data/` | Sheet | Column |
| --- | --- | --- | --- |
| `STRING.py` | `Supplementary File 1.xlsx` | `Table S5` | `Ensembl_ID` |
| `IntAct.py` | `OMIM_gene_mapped_comparison.xlsx` | `Summary` | `Ensembl_ID` |
| `WP.py` | `OMIM_gene_mapped_comparison.xlsx` | `Summary` | `Ensembl_ID` |

Edit `INPUT_FILE`, `SHEET_NAME`, and `ENSEMBL_COLUMN` at the top of each script to match your workbook. To analyze the same gene set across all three resources, point all scripts to the same input.


## Interaction filters

### STRING

`STRING.py` uses human STRING version 12.0 data. Ensembl gene IDs are mapped to STRING protein IDs through the aliases file. An interaction is retained when at least one endpoint maps to an input gene.

Default filters:

- Final combined score: **990–1000**.
- Experimental evidence: `experimental > 0` **or** `experimental_transferred > 0`.

### IntAct

`IntAct.py` maps input Ensembl gene IDs to UniProt accessions and reads the human IntAct MITAB archive. It retains human–human records with usable UniProt identifiers, applies protein-type and detection-method filters, and requires at least one endpoint to match the input set.

The final MI-score range is **0.60–1.00**. The detection-method filter excludes missing methods and records matching specified inferred or predicted labels. This is an exclusion-based filter rather than an exhaustive classification of experimental methods.

### WikiPathways

`WP.py` downloads a human GPML archive from the current WikiPathways release directory. It retains `GeneProduct` and `Protein` data nodes with explicit Ensembl cross-references and requires at least one interaction endpoint to match an input gene.

The accepted interaction annotations are:

- `mim-binding`
- `mim-complex`
- `mim-catalysis`
- `mim-stimulation`
- `mim-inhibition`
- `mim-necessary-stimulation`
- `mim-modification`

No numerical confidence threshold is applied. These annotations describe pathway relationships and do not all imply direct physical binding.



## Outputs

| File | Contents |
| --- | --- |
| `string_ppi_edges.csv` | STRING protein pairs, endpoint identifiers and labels, combined scores, and evidence-channel scores. |
| `string_edge_counts_by_score.csv` | STRING edge counts at each tested combined-score threshold. |
| `unmapped_ensembl_ids_string.csv` | Input genes without STRING mappings, when reported. |
| `intact_ppi_edges.csv` | IntAct protein pairs, endpoint identifiers and labels, MI scores, evidence counts, and publication metadata. |
| `intact_edge_counts_by_score.csv` | IntAct edge counts at each tested MI-score threshold. |
| `unmapped_ensembl_ids.csv` | Input genes without UniProt mappings, when reported. |
| `wikipathways_ppi_edges.csv` | WikiPathways gene pairs, interaction annotations, pathway identifiers and names, and evidence counts. |
| `genes_not_in_wikipathways.csv` | Input genes not encountered among the parsed interaction endpoints, when reported. |


## Downstream analysis
Part 2: Analysis
The data analysis on the genes, interactions is analysis.ipynb
2.1 Check the input genes for sources and cross-reference
2.2 GO enrichment and pathway overpresentation analysis for the input genes
2.3 Network edge sources distribution 
2.4 HotNet2 output data visualization
