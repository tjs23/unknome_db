import sys, os, re, sqlite3, requests, tarfile, gzip

from collections import defaultdict
from urllib import request
from requests.adapters import HTTPAdapter, Retry

from constants import CLUST_ID_CODE, SPECIES_TIEBREAK_PRIORITY, REF_SPECIES, GO_EVIDENCE_WEIGHTS
from constants import SUBCELL_LOC_TERMS, PROTEIN_EXISTENCES, ALT_UNIPROT_ID_URL, SQL_CHUNK_SIZE, OBSOLETE_TAXID
from constants import ORGANISM_DBS, ORGANISM_DB_XREF, UNIPROT_CC_FIELDS, UNIPROT_FT_FIELDS, UNIPROT_SEARCH_URL

import sql_tables
 
# Populate the bulk of data in the Unknome DB using the UniProt website's RESTful API

def _chunked_execute(connection, data_rows, sql_smt, label='rows'):
    """
    If there is a problem you'll know roughly where...
    """    
    cursor = connection.cursor()

    n = len(data_rows)
    
    for i in range(0, n, SQL_CHUNK_SIZE):
        j = min(i+SQL_CHUNK_SIZE, n)
        print(f' .. {label} {i:,} - {j:,}', end='\r')
        cursor.executemany(sql_smt, data_rows[i:j])
    
    print(f' .. {label} {n:,}')
    cursor.close()
    connection.commit()


def _get_uniprot_subcell_locations(terms):

    locations = set()
    
    for term in terms:
        term = term.lower()
        
        for key in SUBCELL_LOC_TERMS:
            if key in term:
                locations.add(SUBCELL_LOC_TERMS[key])
    
    return sorted(locations)


def download_data_files(data_dir, eco_obo_url, panther_class_url_dir, uniprot_go_gaf_url, ncbi_taxonomy_url, buffer_size=2 ** 22):
    
    print('Cached, out-of-date data files may be deleted prior to running this.')
    
    def _buffered_download(url, save_path):
        print(f'Downloading {url} to {save_path}')
        req = request.Request(url)
        
        with request.urlopen(req) as stream_in, open(save_path, 'wb', buffer_size) as file_out:
            response = stream_in.read(buffer_size)
            i = 0
            while response:
                file_out.write(response)
                response = stream_in.read(buffer_size)
                print(f' .. {i}', end='\r')
                i += 1
                
        print('.. done')    
    
    eco_obo_file = os.path.join(data_dir, os.path.basename(eco_obo_url))
    go_gaf_file = os.path.join(data_dir, os.path.basename(uniprot_go_gaf_url))
    
    if os.path.exists(eco_obo_file):
        print(f'Found eco.obo file {eco_obo_file}')
    
    else:
        _buffered_download(eco_obo_url, eco_obo_file)
    
    if os.path.exists(go_gaf_file):
        print(f'Found UniProt GO gaf file {go_gaf_file}')
    
    else:
        _buffered_download(uniprot_go_gaf_url, go_gaf_file)

    with request.urlopen(request.Request(panther_class_url_dir)) as stream_in:
        lines = stream_in.read().decode('utf-8').split('\n')
        
        for line in lines:
            match = re.search('a\s+href="PANTHER(\d.+)_HMM_classifications"', line)
            
            if match:
                panther_version = match.group(1)
                file_name = f'PANTHER{panther_version}_HMM_classifications'
                panther_class_url = os.path.join(panther_class_url_dir, file_name)
        
        print(f'Found PANTHER classifications URL {panther_class_url}')
    
    taxonomy_dir = os.path.join(data_dir, 'taxonomy')
    if not os.path.exists(taxonomy_dir):
        os.mkdir(taxonomy_dir)
    
    tax_names_path = os.path.join(taxonomy_dir, 'names.dmp')
    tax_nodes_path = os.path.join(taxonomy_dir, 'nodes.dmp')
    
    from subprocess import call
    
    if not (os.path.exists(tax_names_path) and os.path.exists(tax_nodes_path)):
        taxonomy_tar_Z = os.path.join(taxonomy_dir, os.path.basename(ncbi_taxonomy_url))
        taxonomy_tar = os.path.splitext(taxonomy_tar_Z)[0]
        _buffered_download(ncbi_taxonomy_url, taxonomy_tar_Z)
        
        print(f'Decompressing  {taxonomy_tar_Z} to {taxonomy_tar}')
        call(['gunzip', '-f', taxonomy_tar_Z])
               
        tar_fileobj = tarfile.open(taxonomy_tar, 'r')
        tar_fileobj.extract('names.dmp', taxonomy_dir)
        tar_fileobj.extract('nodes.dmp', taxonomy_dir)
        
        print(f'Extracted {tax_names_path} and {tax_nodes_path}')
        os.unlink(taxonomy_tar)
     
    return eco_obo_file, go_gaf_file, panther_class_url, f'Panther{panther_version}', tax_names_path, tax_nodes_path
 

def get_panther_name_dict(panther_class_url, cache_file='panther_fam_names.txt'):
     
    if os.path.exists(cache_file):
        print(f'Reading cache {cache_file}')
        
    else:
        n = 0
        with open(cache_file, 'w') as out_file_obj:
            write = out_file_obj.write
        
            req = request.Request(panther_class_url)
 
            with request.urlopen(req) as f:
                response = f.read()
 
                lines = response.decode('utf-8').split('\n')
 
                for line in lines[1:]:
                    line = line.rstrip('\n')
 
                    if not line:
                        continue
 
                    fam_id, fam_name, *null = line.split('\t')
                    write(f'{fam_id}\t{fam_name}\n')
                    n += 1
        
        print(f'Wrote {n:,} lines to {cache_file}')
        
    panther_name_dict = {}
    with open(cache_file) as file_obj:
        for line in file_obj:
            fam_id, fam_name = line.strip().split('\t')
            panther_name_dict[fam_id] = fam_name
            
    return panther_name_dict
    

def get_uniprot_primary_ids(alt_id_url=ALT_UNIPROT_ID_URL, cache_file='uniprot_sec_acc.txt'):
    """
    Helpful for user queries on alt accessions
    """

    print(f'Fetching UniProt secondary accessions')
    
    alt_id_dict = defaultdict(list)
    
    if not os.path.exists(cache_file):
        req = request.Request(alt_id_url)
 
        with request.urlopen(req) as f:
             response = f.read().decode('utf-8')
             
             print(f' .. writing {cache_file}')
             with open(cache_file, 'w') as file_obj:
                 file_obj.write(response)
    
    with open(cache_file) as file_obj:
        lines = file_obj.readlines()
 
    for i, line in enumerate(lines):
        if line.startswith('Secondary AC'):
            break
    
    for line in lines[i+1:]:
        if line.strip():
            sec_id, prim_id = line.split()
            alt_id_dict[prim_id].append(sec_id)
    
    return alt_id_dict
    

def add_ncbi_taxa(db_file_path, names_file, nodes_file):

    connection = sqlite3.connect(db_file_path)

    print('Loading taxonomy names')
    name_dict = {}
    with open(names_file, 'r') as file_obj:
        for line in file_obj:
            tax_id, name, uniq_name, name_class, comment = line.split('|')
            if name_class.strip() == 'scientific name':
                name_dict[int(tax_id)] = name.strip()
    
    smt = 'INSERT INTO NcbiTaxon (tax_id, parent_id, rank, sci_name, embl_code, division,'
    smt += ' inherits_div, genetic_code, inherits_gc, mito_genetic_code, inherits_mito_gc, genbank_hidden, no_seq_data, comments)'
    smt += ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'
    
    cursor = connection.cursor()
    print('Loading taxonomy nodes')
    
    tax_data_rows = []
    
    with open(nodes_file, 'r') as file_obj:
        for line in file_obj:
            data              = line.split('|')
            tax_id            = int(data[0])     # node id in GenBank taxonomy database
            parent_id         = int(data[1])     # parent node id in GenBank taxonomy database
            rank              = data[2].strip()  # rank of this node (superkingdom, kingdom, ...)
            embl_code         = data[3].strip()  # locus-name prefix; not unique
            division          = int(data[4])     # see division.dmp file
            inherits_div      = bool(data[5])    # 1 if node inherits division from parent
            genetic_code      = int(data[6])     # see gencode.dmp file
            inherits_gc       = bool(data[7])    # 1 if node inherits genetic code from parent
            mito_genetic_code = int(data[8])     # see gencode.dmp file
            inherits_mito_gc  = bool(data[9])    # 1 if node inherits mitochondrial gencode from parent
            genbank_hidden    = bool(data[10])   # 1 if name is suppressed in GenBank entry lineage
            no_seq_data       = bool(data[11])   # 1 if this subtree has no sequence data yet
            comments          = data[12].strip()
            sci_name = name_dict.get(tax_id)     # Scientific name from names.dmp file

            row_data = (tax_id, parent_id, rank, sci_name, embl_code, division,
                        inherits_div, genetic_code, inherits_gc, mito_genetic_code,
                        inherits_mito_gc, genbank_hidden, no_seq_data, comments)
            
            tax_data_rows.append(row_data)
            
    cursor.executemany(smt, tax_data_rows)
    cursor.close()
    connection.commit()
    connection.close()
    
    print(f'Inserted {len(tax_data_rows):,} NCBI taxa')
    
    
def add_evidence_codes(db_file_path, eco_file):
    
    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()

    n = 0
    smt = 'INSERT INTO EvidenceCode (eco_id, go_ev_code, name, synonym, knownness_weight) VALUES (?, ?, ?, ?, ?)'
 
    with open(eco_file) as file_obj:
        go_ev_code = None
     
        for line in file_obj:
            if line[:6] == '[Term]':
                go_ev_code = None
            
            if line[:3] == 'id:':
                eco_id = line.split()[1]

            if line[:5] == 'name:':
                name = line[6:-1]
            
            if (line[:5] == 'xref:') and ('GOECO:' in line[6:]):
                go_ev_code = line.split('GOECO:')[1].split()[0]
                synonym = line.split('"')[1]
                
            if go_ev_code and not line.strip(): # end of entry
                knownness_weight = GO_EVIDENCE_WEIGHTS.get(go_ev_code, 0.0)
                cursor.execute(smt, (eco_id, go_ev_code, name, synonym, knownness_weight))
                go_ev_code = None
                n += 1
     
    cursor.close()
    connection.commit()
    connection.close()
                
    print('Inserted %d GO evidence codes' % n)    


def calc_knownness(db_file_path):

    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
    
    print('Calculating knownness and taxonomic lineages')        
    
    evidence_weights = {}
    smt = 'SELECT go_ev_code, knownness_weight FROM EvidenceCode'
    for go_ev_code, knownness_weight in cursor.execute(smt):
        evidence_weights[go_ev_code] = knownness_weight
    
    prot_scores = defaultdict(float)
    prot_go_counts = defaultdict(int)
    smt = 'SELECT go_id, domain, go_ev_code, prot_id FROM ProteinGoTerm'
    for go_id, domain, go_ev_code, prot_id in cursor.execute(smt):
        if domain != 'C':
            prot_scores[prot_id] += evidence_weights[go_ev_code]
     
        prot_go_counts[prot_id] += 1

    prot_tax_id = {}
    prot_names = {}
    smt = 'SELECT id, name, panther_db_name, species_id FROM Protein'
    for prot_id, name, panther_db_name, tax_id in cursor.execute(smt):
        prot_tax_id[prot_id] = tax_id
        prot_names[prot_id] = panther_db_name or name
        
    print(f'Insert protein knownness')
    
    prot_known_rows = []
    smt = 'UPDATE Protein SET knownness=? WHERE id=?'
    for prot_id, knownness in prot_scores.items():
        prot_known_rows.append((knownness, prot_id))
        
        if len(prot_known_rows)> SQL_CHUNK_SIZE:
            cursor.executemany(smt, prot_known_rows)
            prot_known_rows = []
    
    if prot_known_rows:
        cursor.executemany(smt, prot_known_rows)
    
    print(f'Calculate cluster knownness')
        
    cluster_scores = defaultdict(float)
    cluster_scores_best = defaultdict(str)
    cluster_counts = defaultdict(float)
    cluster_counts_best = defaultdict(str)
    cluster_lineages = defaultdict(set)

    smt = 'SELECT clust_id, prot_id FROM OrthoClusterProteins'
    for clust_id, prot_id in cursor.execute(smt):
        score = prot_scores[prot_id]
        tax_id = prot_tax_id[prot_id]
        
        if tax_id in REF_SPECIES:
            if score > cluster_scores[clust_id]:
                cluster_scores[clust_id] = score
                cluster_scores_best[clust_id] = prot_id
 
            go_count = prot_go_counts[prot_id] + SPECIES_TIEBREAK_PRIORITY.get(tax_id, 0.0)
            if go_count > cluster_counts[clust_id]:
                cluster_counts[clust_id] = go_count
                cluster_counts_best[clust_id] = prot_id
                
        cluster_lineages[clust_id].add(tax_id)
        del prot_scores[prot_id]
    
    ortho_cluster_rows = []
    print(f'Insert cluster knownness')
    
    smt = 'UPDATE OrthoCluster SET name=?, knownness=?, best_known=?, tax_lineages=? WHERE id=?'
    for clust_id, knownness in cluster_scores.items():
        best_known = cluster_scores_best[clust_id] or cluster_counts_best[clust_id]
        
        if best_known:
            name = prot_names[best_known]
        else:
            name = 'Unknown'
            
        tax_lineages = len(cluster_lineages[clust_id])
        row_data = (name, knownness, best_known, tax_lineages, clust_id)
        ortho_cluster_rows.append(row_data)
        
        if len(ortho_cluster_rows) > SQL_CHUNK_SIZE:
            cursor.executemany(smt, ortho_cluster_rows)
            ortho_cluster_rows = []            
    
    if ortho_cluster_rows:
        cursor.executemany(smt, ortho_cluster_rows)    

    cursor.close()
    connection.commit()
    connection.close()

    print(' .. done')        
    
        
def add_uniprot_proteins_panther_clusters(panther_class_url, odb_source='Panther17', verbose=False):
    
    panther_name_dict = get_panther_name_dict(panther_class_url)

    print('Adding proteins and orthologue groups')
    
    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
    
    # From UniProt examples

    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_uniprot_next_link(headers):
     if "Link" in headers:
         match = re_next_link.match(headers["Link"])
         if match:
             return match.group(1)

    def get_uniprot_batch(batch_url):
        while batch_url:
             response = session.get(batch_url)
             response.raise_for_status()
             yield response
             batch_url = get_uniprot_next_link(response.headers)
   
    n_organism_dbs = len(ORGANISM_DBS)
    
    main_fields = ('id','accession','protein_name', 'gene_primary', 'gene_synonym','gene_oln', 'gene_orf', 'organism_name', 'organism_id',
                   'lineage','xref_embl', 'xref_pfam', 'xref_interpro','xref_panther', 'go', 'reviewed', 'protein_existence', 'sequence_version', 'sequence')
    
    fields = main_fields + tuple([ORGANISM_DB_XREF[odb] for odb in ORGANISM_DBS]) + UNIPROT_CC_FIELDS + UNIPROT_FT_FIELDS
    fields = ','.join(fields)
   
    uniprot_url = f"{UNIPROT_SEARCH_URL}query=database:panther&format=tsv&size=500&fields={fields}"
     
    alt_acc_id_dict = get_uniprot_primary_ids()
    
    p1 = len(main_fields)
    p2 = p1+n_organism_dbs
    p3 = p2+len(UNIPROT_CC_FIELDS)
    
    tax_info_dict = {}
    go_term_dict = {}
    cluster_ids = {}
    
    duplic_acc = {}
    prot_accession_rows = []
    prot_rows = []
    prot_cluster_rows = []
    cluster_rows = []
    
    progress = 0
    
    clust_add_smt = 'INSERT INTO OrthoCluster (id, db_source, db_id, knownness, name, clust_group) VALUES (?, ?, ?, ?, ?, ?)'
    prot_add_smt  = 'INSERT INTO Protein (id, accessions, name, gene_name, orf_name, species, species_id, species_aka,'
    prot_add_smt += 'tax_id, tax_name, tax_lineage, subcell_locs, pubmed_ids, embl_ids, pfam_ids, interpro_ids, '
    prot_add_smt += 'organism_db_xref, organism_db, panther_db_xref, panther_db_name, clust_id, features, fasta_header, sequence)'
    prot_add_smt += ' VALUES (%s)' % ', '.join(['?'] * (prot_add_smt.count(',')+1))
    prot_acc_smt  = 'INSERT INTO ProteinAccession (accession, prot_id) VALUES (?, ?)'
    prot_clust_smt = 'INSERT INTO OrthoClusterProteins (clust_id, prot_id) VALUES (?, ?)'

    for batch in get_uniprot_batch(uniprot_url):
        lines = batch.text.splitlines()
        
        if not progress:
            print(lines[0])
        
        for line in lines[1:]:
            row = line.split('\t')
            prot_id, accession, protein_name, gene_primary, gene_synonym, gene_oln, gene_orf, organism_name, tax_id = row[0:9]
            lineage, xref_embl, xref_pfam, xref_interpro, xref_panther, go, reviewed, protein_existence, sequence_version, sequence = row[9:p1]
            org_db_refs = row[p1:p2]
            cc_vals = row[p2:p3]
            ft_vals = row[p3:]
            
            tax_id = int(tax_id)
            if tax_id in OBSOLETE_TAXID:
              tax_id = OBSOLETE_TAXID[tax_id]
            
            gene_name = gene_primary or gene_oln or gene_orf 
            pe = PROTEIN_EXISTENCES[protein_existence]
            db_section = 'sp' if (reviewed == 'reviewed') else 'tr'
            
            fasta_header = f'{db_section}|{accession}|{prot_id} {protein_name} OS={organism_name} OX={tax_id} GN={gene_name} PE={pe} SV={sequence_version}'
            
            xref_embl = xref_embl.strip(';') or None
            xref_pfam = xref_pfam.strip(';') or None
            xref_interpro = xref_interpro.strip(';') or None
            #xref_orthodb = xref_orthodb.strip(';') or None            
            #xref_treefam = xref_treefam.strip(';') or None            
            #xref_eggnog = xref_eggnog.strip(';') or None
            #xref_ko = xref_ko.strip(';') or None

            xref_panther = xref_panther.strip(';')
            
            if ';' in xref_panther:
                panther1, panther2 = xref_panther.split(';')
 
                if ':SF' in panther1:
                    panther_fam, panther_sf = panther2, panther1
                else:
                    panther_fam, panther_sf = panther1, panther2
            
            else:
                panther_fam = xref_panther
            
            scl = cc_vals[6]
            scl = scl.strip('.')
            
            if scl:
                scl.replace('SUBCELLULAR LOCATION: ', '')
            
                if 'Note=' in scl:
                    scl = scl.split('Note=')[0].strip()
                        
                loc_terms = [l.strip().split(';')[0].split(' {ECO:')[0] for l in scl.split('. ')]
                loc_terms = [l.split(', ')[0].split('.')[0].split(':')[-1].strip() for l in loc_terms if l]

                subcell_locs = ';'.join(_get_uniprot_subcell_locations(loc_terms))
            
            else:
                subcell_locs = None
            
            pubmed_ids = []
            for cc in cc_vals:
                for pmid in re.findall('PubMed:\d+', cc):
                    pubmed_ids.append(pmid.split(':')[1])
            
            pubmed_ids = ';'.join(pubmed_ids)
            
            for odb, xref in zip(ORGANISM_DBS, org_db_refs):
                if xref:
                    xref = ';'.join([x.split()[0] for x in xref.strip(';').split(';')])
                    organism_db_xref = xref
                    organism_db = odb
                    break
            
            # Get taxon ID
            if tax_id in tax_info_dict: # Cached previously
                species, species_id, species_aka, tax_name = tax_info_dict[tax_id]
            
            else:
                smt = 'SELECT parent_id, rank, sci_name FROM NcbiTaxon WHERE tax_id=?'
                row =  cursor.execute(smt, [tax_id]).fetchone()
                
                if row:
                    tax_parent_id, tax_rank, tax_name =row
 
                    # Get parent species
                    species_id = tax_id
                    while tax_rank != 'species':
                        species_id = tax_parent_id
                        tax_parent_id, tax_rank, null = cursor.execute(smt, [tax_parent_id]).fetchone()
 
                    species = ' '.join(organism_name.split()[:2])
                    if '(' in species:
                        species_aka = organism_name.split('(')[1].strip()[:-1]
                    else:
                        species_aka = None
 
                    if verbose:
                        print(f'Found species {tax_id}/{species_id} {species}/{species_aka} {tax_name}')
                    tax_info_dict[tax_id] = species, species_id, species_aka, tax_name
                
                else:
                    print(f'No species info for tax_id {tax_id} at {prot_id}: NCBI taxonomy .dmp files could be out-of-date if this occurs often')
                    continue
                
            accessions = [accession] + alt_acc_id_dict[accession]
            
            for acc in accessions:
                if acc in    duplic_acc:
                    if verbose:
                        print(f'Repeated protein accession {acc} refers to {prot_id} and {duplic_acc[acc]}')
                    continue
             
                duplic_acc[acc] = prot_id
                prot_accession_rows.append((acc, prot_id))
            
            accessions = accessions = ';'.join(accessions)
            

            features = {}
            for ft_val in ft_vals:
                if ft_val:
                    key = ft_val.split()[0]
                    count = ft_val.count(key)
                    features[key] = count
             
            features = ';'.join(['%s:%d' % (x, features[x]) for x in sorted(features)])
            
            if panther_fam not in panther_name_dict:
                if verbose:
                    print(f'PANTHER family {panther_fam} on entry {prot_id} ({species}) does not exist')
                continue
                
            panther_db_name = panther_name_dict[panther_fam]
            
            go = go.strip()
            
            if go:
                for item in go.strip(']').split(']; '):
                    term, go_id = item.split('[GO:')
                    term = term.strip()
                    go_id = 'GO:' + go_id
                    go_term_dict[go_id] = term
                                
            if panther_fam in cluster_ids:
                clust_id = cluster_ids[panther_fam]
            
            else:
                clust_id = len(cluster_ids)
                cluster_ids[panther_fam] = clust_id
                cluster_rows.append((CLUST_ID_CODE % clust_id, odb_source, panther_fam, None, panther_db_name, 0))
                #(id, db_source, db_id, knownness, name, clust_group) 

            prot_cluster_rows.append((clust_id, prot_id)) # Added after clusters all added
         
            prot_rows.append((prot_id, accessions, protein_name, gene_name, gene_orf, species, species_id,
                              species_aka, tax_id, tax_name, lineage, subcell_locs, pubmed_ids, xref_embl,
                              xref_pfam, xref_interpro, organism_db_xref, organism_db, panther_fam,
                              panther_db_name, clust_id, features, fasta_header, sequence))
            
            if len(prot_rows) > SQL_CHUNK_SIZE:
                if cluster_rows:
                    cursor.executemany(clust_add_smt, cluster_rows)
                    cluster_rows = []
            
                cursor.executemany(prot_add_smt, prot_rows)
                prot_rows = []
                
                cursor.executemany(prot_acc_smt, prot_accession_rows)
                prot_accession_rows = []
                
                cursor.executemany(prot_clust_smt, prot_cluster_rows) # cluster and prot IDs alreadt inserted
                prot_cluster_rows = []
                
                
        progress += len(lines[1:])
        
        print(f'{progress}')
    
    # Any remaining final, part chunks
    
    if cluster_rows:
        cursor.executemany(clust_add_smt, cluster_rows)
                
    if prot_rows:
        cursor.executemany(prot_add_smt, prot_rows)

    if prot_accession_rows:
        cursor.executemany(prot_acc_smt, prot_accession_rows)

    if prot_cluster_rows:
        cursor.executemany(prot_clust_smt, prot_cluster_rows)
    
    cursor.close()
    connection.commit()
    connection.close()

    print(' .. done')        
        
    return go_term_dict
    
    
def add_go_terms_with_dates(db_file_path, uniprot_go_gaf_path, go_term_dict):

    print('Adding protein GO term evidence codes and dates')        

    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
 
    pid_map = {}
    prot_pubmed = {}
    for pid, accessions, pubmed_ids in cursor.execute('SELECT id, accessions, pubmed_ids FROM Protein'):
        if pubmed_ids:
            prot_pubmed[pid] = pubmed_ids.split(';')
        
        for acc in accessions.split(';'):
            pid_map[acc] = pid    
 
    cursor.close()
    
    n = 0
    protein_go_terms = {}
    protein_pubmed_ids = []
    
    with gzip.open(uniprot_go_gaf_path, 'r',    2**12) as file_obj:
        for line in file_obj:
            if line[0] == '!':
                continue
            
            data = line.split('\t')
            go_id = data[4]
            ref = data[5]
            go_ev_code = data[6]
            prot_id = pid_map.get(data[1])
            
            if prot_id: # In the PATHER set
                date_str = data[13]
                year = int(date_str[:4])
                month = int(date_str[4:6])
                day = int(date_str[6:])
                source = data[7].split(':')[0]
                go_domain = data[8] # F/P/C
                term = go_term_dict[go_id]
                
                key = (go_id, go_ev_code, prot_id)
                
                if key in protein_go_terms:
                    year_prev, month_prev, day_prev, src_prev, dom_prev, term_prev = protein_go_terms[key]
                    
                    # Swap only if newer
                    if (year > year_prev) or ((year_prev == year) and ((month > month_prev) or ((month_prev == month) and (day > day_prev)))):
                        protein_go_terms[key] = (year, month, day, source, go_domain, term)                
                                            
                else:
                    protein_go_terms[key] = (year, month, day, source, go_domain, term)                
                
                if ref[:5] == 'PMID:':
                    pmids = prot_pubmed.get(prot_id, [])
                    pmid = ref[5:]
                    
                    if pmid not in pmids:
                        pmids.append(pmid)

                        pmids = ';'.join(sorted(set(pmids)))
                        protein_pubmed_ids.append((pmids, prot_id))
                
                n += 1

    protein_go_terms = [k+v for k, v in protein_go_terms.items()]
    
    print(f"    ..found {n:,}")
    
    print(f'Adding {len(protein_go_terms):,} protein GO terms')        
    smt = 'INSERT INTO ProteinGoTerm (go_id, go_ev_code, prot_id, date_year, date_month, date_day, project, domain, term) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)'
    _chunked_execute(connection, protein_go_terms, smt, 'GO terms')
    
    print(f'Updating {len(protein_pubmed_ids):,} protein PubMed IDs')        
    smt = 'UPDATE Protein SET pubmed_ids=? WHERE id=?'
    _chunked_execute(connection, protein_pubmed_ids, smt, 'PubMed refs')
    
    connection.close()
    
    print(' .. done')        
       
    
def make_indices(db_file_path):

    print('Adding DB indices')        

    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
    
    # OrthoCluster
    cursor.execute("CREATE INDEX ortho_cluster_knownness ON OrthoCluster(knownness)") 
    cursor.execute("CREATE INDEX ortho_cluster_lineages ON OrthoCluster(tax_lineages)") 
    cursor.execute("CREATE INDEX ortho_cluster_group ON OrthoCluster(clust_group)") 
    
    # Protein
    cursor.execute("CREATE INDEX protein_species ON Protein(species)") 
    cursor.execute("CREATE INDEX protein_name ON Protein(name)") 
    cursor.execute("CREATE INDEX protein_gene_name ON Protein(gene_name)") 
    cursor.execute("CREATE INDEX protein_orf_name ON Protein(orf_name)") 
    cursor.execute("CREATE INDEX protein_clust_id ON Protein(clust_id)") 
    
    # OrthoClusterProteins
    cursor.execute("CREATE INDEX cluster_proteins_clust_id ON OrthoClusterProteins(clust_id)") 
    
    # ProteinGoTerm
    cursor.execute("CREATE INDEX protein_go_term_prot_id ON ProteinGoTerm(prot_id)") 
    cursor.execute("CREATE INDEX protein_go_term_go_id ON ProteinGoTerm(go_id)") 

    cursor.close()
    connection.commit()
    connection.close()
    
    print(' .. done')        
    

def make_tables(db_file_path):

    print('Making DB tables')        
    
    if os.path.exists(db_file_path):
        os.unlink(db_file_path)

    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
    
    cursor.execute(sql_tables.TABLE_EVIDENCE_CODE)
    cursor.execute(sql_tables.TABLE_NCBI_TAXON)
    cursor.execute(sql_tables.TABLE_CLUSTER_GROUP)
    cursor.execute(sql_tables.TABLE_CLUSTER)
    cursor.execute(sql_tables.TABLE_PROTEIN)
    cursor.execute(sql_tables.TABLE_CLUSTER_PROTEINS)
    cursor.execute(sql_tables.TABLE_PROTEIN_GO_TERM)
    cursor.execute(sql_tables.TABLE_PROTEIN_ACCESSION)
    
    cursor.close()
    connection.commit()
    
    print(' .. done')        
    

if __name__ == '__main__':
    
    # Directory for storage of data used to fill DB    
    data_path = '/data/unknome/'
    
    # File path to the UNknome DB in sqlite3 format
    db_file_path = '/data/unknome/unknome_db_Sep_2023.sqlite'
    
    # URL for Gene Ontology eveidence codes
    eco_obo_url = 'https://github.com/evidenceontology/evidenceontology/blob/master/eco.obo'
    
    # URL for the DIRECTORY containing the latest PANTHER HMM classifications file
    panther_class_url_dir = "http://data.pantherdb.org/ftp/hmm_classifications/current_release/"
    
    # URL for GeneOntology to UniProt GOA .gaf file. 
    uniprot_go_gaf_url = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
    
    # URL for TAR GZIPPEd archive containing NCBI taxonomy .dmp files
    ncbi_taxonomy_url = 'https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z'
    
    # Download data files to data dir where they do not exist
    paths = download_data_files(data_path, eco_obo_url, panther_class_url_dir, uniprot_go_gaf_url, ncbi_taxonomy_url)
    eco_file_path, uniprot_go_gaf_path, panther_class_url, panther_version, tax_names_path, tax_nodes_path = paths
    
    # Make a fresh databse file and create tables
    make_tables(db_file_path)
    
    # Insert GO Evidence Codes from eco.obo
    add_evidence_codes(db_file_path, eco_file_path)
    
    # Insert NCBI taxonomy data; handy for working out larger clade membersships, and esp. which species a strain belongs to
    add_ncbi_taxa(db_file_path, tax_names_path, tax_nodes_path)
    
    # Insert UniProt protein infor where there are PANTHER links; fills protein and orthologue cluster tables
    go_term_dict = add_uniprot_proteins_panther_clusters(panther_class_url, panther_version)
    
    # Insert dates for protein GO terms using .gaf file
    add_go_terms_with_dates(db_file_path, uniprot_go_gaf_path, go_term_dict)             
 
    # Calculate knownness values for clusters and proteins given protein GO links
    calc_knownness(db_file_path)
    
    # Make DB indices for faster lookups
    make_indices(db_file_path)             

