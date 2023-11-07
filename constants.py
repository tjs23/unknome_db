
FILE_END = 'END\n'

CLUST_ID_CODE = 'UKP%05d'

# K12 Ecoli, Human, mouse,yeast,pombe, arabidopsis,fly,worm,fish,frog
SPECIES_TIEBREAK_PRIORITY = {562:0.95, 9696:0.9, 10090:0.8, 4932:0.7, 4896:0.6, 3702:0.5, 7227:0.4,  6239:0.3, 7955:0.2, 8364:0.1}
 
# http://geneontology.org/page/reference-genome-annotation-project 
REF_SPECIES = {3702, 6239, 7955, 44689, 7227, 562, 9031, 9606, 10090, 10116, 4932, 4896}

GO_EVIDENCE_WEIGHTS = {"EXP": 0.8, "IDA": 0.8, "IPI": 0.8, "IMP": 0.8,
                       "IGI": 0.8, "IEP": 0.8, "ISS": 0.5, "ISO": 0.5,
                       "ISA": 0.5, "ISM": 0.5, "IGC": 0.3, "RCA": 0.6,
                       "TAS": 0.9, "NAS": 0.6,  "IC": 1.0,  "ND": 0.0,
                       "IEA": 0.0,  "NR": 0.0, "IRD": 0.0, "IKR": 0.0,
                       "IBA": 0.5, "IBD": 0.5, "HDA": 0.5, "HMP": 0.5,
                       "HEP": 0.5, "HGI": 0.5}

SUBCELL_LOC_TERMS = {'vesicle':'ves', 'endoplasmic':'er',
                     'sarcoplasmic':'er', 'microsome':'er',
                     'endosome':'end', 'golgi':'gol',
                     'lysosome':'lys', 'mitochondri' :'mit',
                     'chromosome':'nuc', 'nucleus':'nuc',
                     'plastid':'pla', 'thylakoid':'pla',
                     'cytoplasm':'cyt', 'droplet':'cyt',
                     'extracellular':'ext', 'periplasm':'ext',
                     'fimbrium':'ext', 'secreted':'ext', 'virion':'ext',
                     'cell membrane':'pm', 'flagellum':'pm',
                     'cell inner membrane':'pm', 'cell outer membrane':'pm',
                     'cell envelope':'pm', 'cell junction':'pm',
                     'cell surface':'pm', 
                     'cleavage furrow':'pm', 'midbody':'pm',
                     'glycosome':'oth', 'glyoxysome':'oth',
                     'peroxisome':'oth', 'preautophagosomal':'oth',
                     'melanosome':'oth', 'prospore':'oth',
                     'spore':'oth', 'vacuole':'oth', 'vacuolar':'oth'}

PROTEIN_EXISTENCES = {'Evidence at protein level':1,
                      'Evidence at transcript level':2,
                      'Inferred from homology':3,
                      'Predicted':4,
                      'Uncertain':5}

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search?"

ALT_UNIPROT_ID_URL = 'https://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/complete/docs/sec%5Fac.txt'

SQL_CHUNK_SIZE = 4000

OBSOLETE_TAXID = {342610:3042615, # c.f. UniProt:CCA_PSEA6
                  168810:1718,    # c.f. UniProt:ARGJ_CORCT
                  } 

ORGANISM_DBS = ['FlyBase', 'MaizeGDB', 'RGD', 'SGD', 'TAIR', 'PomBase', 'WormBase', 'Xenbase', 'ZFIN', 'dictyBase',
                'Ensembl', 'EnsemblBacteria', 'EnsemblFungi', 'EnsemblMetazoa', 'EnsemblPlants', 'EnsemblProtists',
                'AG', 'ArachnoServer', 'Araport', 'CGD', 'CTD', 'ConoServer', 'DisGeNET', 'EchoBASE', 'GeneCards',
                'GeneReviews', 'HGNC', 'HPA', 'LegioList', 'Leproma', 'MGI', 'MIM', 'MalaCards', 'NIAGADS', 'OpenTargets',
                'Orphanet', 'PharmGKB', 'PseudoCAP', 'TubercuList', 'VEuPathDB', 'VGNC', 'euHCVdb', 'neXtProt', ] # Priority sorted 

ORGANISM_DB_XREF = {'Ensembl':'xref_ensembl', 'EnsemblBacteria':'xref_ensemblbacteria', 'EnsemblFungi':'xref_ensemblfungi',
                    'EnsemblMetazoa':'xref_ensemblmetazoa', 'EnsemblPlants':'xref_ensemblplants', 'EnsemblProtists':'xref_ensemblprotists',
                    'AG':'xref_agr', 'ArachnoServer':'xref_arachnoserver', 'Araport':'xref_araport', 'CGD':'xref_cgd', 'CTD':'xref_ctd',
                    'ConoServer':'xref_conoserver', 'DisGeNET':'xref_disgenet', 'EchoBASE':'xref_echobase', 'FlyBase':'xref_flybase',
                    'GeneCards':'xref_genecards', 'GeneReviews':'xref_genereviews', 'HGNC':'xref_hgnc', 'HPA':'xref_hpa',
                    'LegioList':'xref_legiolist', 'Leproma':'xref_leproma', 'MGI':'xref_mgi', 'MIM':'xref_mim', 'MaizeGDB':'xref_maizegdb',
                    'MalaCards':'xref_malacards', 'NIAGADS':'xref_niagads', 'OpenTargets':'xref_opentargets', 'Orphanet':'xref_orphanet',
                    'PharmGKB':'xref_pharmgkb', 'PomBase':'xref_pombase', 'PseudoCAP':'xref_pseudocap', 'RGD':'xref_rgd', 'SGD':'xref_sgd',
                    'TAIR':'xref_tair', 'TubercuList':'xref_tuberculist', 'VEuPathDB':'xref_veupathdb', 'VGNC':'xref_vgnc',
                    'WormBase':'xref_wormbase', 'Xenbase':'xref_xenbase', 'ZFIN':'xref_zfin', 'dictyBase':'xref_dictybase',
                    'euHCVdb':'xref_euhcvdb', 'neXtProt':'xref_nextprot'}

UNIPROT_CC_FIELDS = ('cc_function', 'cc_catalytic_activity', 'cc_cofactor', 'cc_activity_regulation','cc_subunit', 
                     'cc_pathway', 'cc_subcellular_location', 'cc_tissue_specificity', 'cc_developmental_stage',
                     'cc_induction', 'cc_domain', 'cc_ptm', 'cc_rna_editing', 'cc_mass_spectrometry', 'cc_polymorphism',
                     'cc_disease', 'cc_disruption_phenotype',
                     'cc_allergen', 'cc_toxic_dose', 'cc_biotechnology', 'cc_pharmaceutical',
                     'cc_miscellaneous', 'cc_similarity', 'cc_caution', 'cc_sequence_caution')

UNIPROT_FT_FIELDS = ('ft_act_site','ft_binding','ft_carbohyd','ft_chain','ft_coiled','ft_compbias',
                     'ft_conflict','ft_crosslnk','ft_disulfid','ft_dna_bind','ft_domain','ft_helix',
                     'ft_init_met','ft_intramem','ft_lipid','ft_mod_res','ft_motif','ft_mutagen',
                     'ft_non_cons','ft_non_std','ft_non_ter','ft_peptide','ft_propep','ft_region',
                     'ft_repeat','ft_signal','ft_site','ft_topo_dom','ft_transit','ft_transmem','ft_turn',
                     'ft_unsure','ft_variant','ft_var_seq','ft_zn_fing')
