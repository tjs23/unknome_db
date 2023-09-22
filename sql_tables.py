TABLE_NCBI_TAXON = """
CREATE TABLE NcbiTaxon (
    tax_id INT,
    parent_id INT,
    rank VARCHAR(16) NOT NULL,
    sci_name TEXT,
    embl_code VARCHAR(2),
    division INT,
    inherits_div BOOLEAN,
    genetic_code INT,
    inherits_gc BOOLEAN,
    mito_genetic_code INT,
    inherits_mito_gc BOOLEAN,
    genbank_hidden BOOLEAN,
    no_seq_data BOOLEAN,
    comments TEXT,
    PRIMARY KEY (tax_id)
);"""

TABLE_EVIDENCE_CODE = """
CREATE TABLE EvidenceCode (
    go_ev_code VARCHAR(4),
    eco_id VARCHAR(12),
    name TEXT,
    synonym TEXT,
    knownness_weight FLOAT,
    PRIMARY KEY (go_ev_code)
);"""

TABLE_CLUSTER_GROUP = """
CREATE TABLE OrthoClusterGroup (
    id INT,
    name VARCHAR(16),
    description TEXT,
    PRIMARY KEY (id)
);"""

TABLE_CLUSTER = """
CREATE TABLE OrthoCluster (
    id VARCHAR(16),
    db_source TEXT NOT NULL,
    db_id TEXT NOT NULL,
    knownness FLOAT,
    name TEXT,
    clust_group INT,
    best_known VARCHAR(16),
    tax_lineages INT,
    FOREIGN KEY (clust_group) REFERENCES ClusterGroup(id),
    PRIMARY KEY (id)
);"""

TABLE_PROTEIN = """
CREATE TABLE Protein (
    id VARCHAR(16),
    accessions TEXT,
    name TEXT,
    gene_name TEXT,
    orf_name TEXT,
    subcell_locs TEXT,
    fasta_header TEXT,
    species TEXT,
    species_aka TEXT,
    species_id INT, 
    tax_id INT, 
    tax_name TEXT, 
    tax_lineage TEXT, 
    pubmed_ids TEXT,
    embl_ids TEXT,
    pfam_ids TEXT,
    interpro_ids TEXT,
    panther_db_xref VARCHAR(16),
    panther_db_name TEXT,
    organism_db_xref VARCHAR(16),
    organism_db_name TEXT,
    organism_db VARCHAR(16),
    clust_id VARCHAR(16),
    knownness FLOAT,
    features TEXT,
    sequence TEXT,
    align_seq TEXT,
    FOREIGN KEY (clust_id) REFERENCES Cluster(id),
    PRIMARY KEY (id)
);"""

TABLE_CLUSTER_PROTEINS = """
CREATE TABLE OrthoClusterProteins (
    id VARCHAR(16),
    clust_id VARCHAR(16),
    prot_id VARCHAR(16),
    FOREIGN KEY (prot_id) REFERENCES Protein(id),
    FOREIGN KEY (clust_id) REFERENCES Cluster(id),
    PRIMARY KEY (id)
);"""

TABLE_PROTEIN_GO_TERM = """
CREATE TABLE ProteinGoTerm (
    go_id VARCHAR(16) NOT NULL,
    domain CHARACTER(1),
    term TEXT,
    go_ev_code VARCHAR(4),
    project VARCHAR(16),
    prot_id VARCHAR(16),
    date_year INT,
    date_month INT,
    date_day INT,
    FOREIGN KEY (prot_id) REFERENCES Protein(id), 
    FOREIGN KEY (go_ev_code) REFERENCES EvidenceCode(go_ev_code), 
    PRIMARY KEY (go_id, go_ev_code, prot_id)
);"""

TABLE_PROTEIN_ACCESSION = """
CREATE TABLE ProteinAccession (
    prot_id VARCHAR(16),
    accession VARCHAR(16),
    FOREIGN KEY (prot_id) REFERENCES Protein(id),
    PRIMARY KEY (accession)
);"""

