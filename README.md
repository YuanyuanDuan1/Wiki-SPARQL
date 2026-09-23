# Genetic epilepsy biomolecular interaction network construction and clustering

Part 1: building the Interaction Network
Input: Ensembl IDs from input file; output: interaction edges from STRING, IntAct, and WikiPathway

1.1 STRING data query 
Variable: combined score(0-1000)
Output columns: 
source_label	source_string_id	source_ensembl	target_string_id	target_label	target_ensembl	combined_score	neighborhood	neighborhood_transferred	fusion	cooccurence	homology	coexpression	coexpression_transferred	experimental	experimental_transferred	database	database_transferred	textmining	textmining_transferred	source_database

1.2 IntAct data query 
Variable: MI score(0-1)
Output columns: source_label	source_uniprot	source_ensembl	target_ensembl	target_label	target_uniprot	intact_miscore	evidence_count	detection_method	interaction_type	publication	publication_count	intact_interaction_id	first_author	source_database

1.3 WikiPathway
variable: None
Only get the MIM interactions from for the input lists and the metabolite/chemical interactions only if any input gene is present in the pathway.
output columns:source_label	source_ensembl	target_label	target_ensembl	interaction_type	pathway_id	pathway_count	pathway_name	evidence_count	source_database


Part 2: Analysis
The data analysis on the genes, interactions is analysis.ipynb
2.1 Check the input genes for sources and cross-reference
2.2 GO enrichment and pathway overpresentation analysis for the input genes
2.3 Network edge sources distribution 
2.4 HotNet2 output data visualization
