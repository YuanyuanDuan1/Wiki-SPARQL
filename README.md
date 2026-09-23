# Genetic epilepsy biomolecular interaction network construction and clustering

Part 1: building the Interaction Network
Input: Ensembl IDs from input file; output: interaction edges from STRING, IntAct, and WikiPathway
1.1 STRING data query 

Output columns: 
source_label	source_string_id	source_ensembl	target_string_id	target_label	target_ensembl	combined_score	neighborhood	neighborhood_transferred	fusion	cooccurence	homology	coexpression	coexpression_transferred	experimental	experimental_transferred	database	database_transferred	textmining	textmining_transferred	source_database

1.2 IntAct data query 

Output columns: source_label	source_uniprot	source_ensembl	target_ensembl	target_label	target_uniprot	intact_miscore	evidence_count	detection_method	interaction_type	publication	publication_count	intact_interaction_id	first_author	source_database

1.3


Part 2: Analyze the network
