import os, sqlite3, subprocess, time, glob
from collections import defaultdict
from plot_summary_charts import *

BASE_DIR = os.path.dirname(__file__)
UNKNOME_DB_FILE_PATH = os.path.join(BASE_DIR, 'unknome_db.sqlite')
TSV_PROT_TABLE_FILE = 'unknome_protein_table_{}.tsv'
TSV_CLUST_TABLE_FILE = 'unknome_cluster_table_{}.tsv'

REF_SPECIES = {3702, 6239, 7955, 44689, 7227, 562, 9031, 9606, 10090, 10116, 4932, 4896}

def _db_date_str(db_file_path):

  real_db_path = os.path.realpath(db_file_path)
  mod_epoch_time = os.path.getmtime(real_db_path)
  db_mod_date = time.strftime('%d_%B_%Y', time.localtime(mod_epoch_time))
 
  return db_mod_date
  
  
def _make_summary_protein_table(download_dir, db_file_path=UNKNOME_DB_FILE_PATH):
  
  file_path = os.path.join(download_dir, TSV_PROT_TABLE_FILE.format(_db_date_str(db_file_path)))

  print(f'Making summary Protein TSV table {file_path}')
  
  connection = sqlite3.connect(db_file_path)
  cursor = connection.cursor()

  smt = 'SELECT id, accessions, knownness, gene_name, name, tax_id, species, clust_id, panther_db_xref FROM Protein'
  head = ['uniprot_name','uniprot_accessions','knownness','gene_name','protein_name','taxon_id','species','cluster_id','panther_group']
  
  with open(file_path, 'w') as file_obj:
    write = file_obj.write
    
    line = '\t'.join(head) + '\n'
    write(line)

    for pid, accessions, knownness, gene_name, name, tax_id, species, clust_id, panther_db_xref in cursor.execute(smt):
      line = '{}\t{}\t{:.3f}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pid, accessions or '', knownness or 0.0, gene_name or '', name or '', tax_id, species, clust_id, panther_db_xref)
      write(line)
 
  cursor.close()
  connection.close()
  
  return file_path
  

def _make_summary_cluster_table(download_dir, db_file_path=UNKNOME_DB_FILE_PATH):
  
  file_path = os.path.join(download_dir, TSV_CLUST_TABLE_FILE.format(_db_date_str(db_file_path)))
  
  print(f'Making summary Cluster TSV table {file_path}')
  
  connection = sqlite3.connect(db_file_path)
  cursor = connection.cursor()
  
  prot_info = defaultdict(list)
  smt = 'SELECT id, gene_name, name, species, species_id, clust_id, organism_db, organism_db_xref FROM Protein'
  
  for pid, gene_name, name, species, species_id, clust_id, organism_db, organism_db_xref in cursor.execute(smt):
    prot_info[clust_id].append((pid, gene_name, name, species_id, species, organism_db, organism_db_xref))
  
  smt = 'SELECT id, db_id, knownness, name, tax_lineages, best_known FROM OrthoCluster ORDER BY knownness ASC'
  
  head = ['cluster_id','panther_id','cluster_name','knownness','num_proteins','num_species',
          'best_known_protein_id','best_known_protein_gene','best_known_protein_name',
          'key_protein_ids','key_protein_xrefs']
  
  with open(file_path, 'w') as file_obj:

    write = file_obj.write
    line = '\t'.join(head) + '\n'
    write(line)

    for cid, db_id, knownness, name, tax_lineages, best_known in cursor.execute(smt):
      
      if not knownness:
        knownness = 0.0
        
      num_prots = len(prot_info[cid])
      best_known_pn = ''
      best_known_gn = ''
      ref_pn = ''
      ref_gn = ''
      ref_id = ''
      cspecies = set()
      key_proteins = defaultdict(list)
      key_protein_ids = set()
      
      for pid, gene_name, prot_name, species_id, species, organism_db, organism_db_xref in prot_info[cid]:
        
        if pid == best_known:
          best_known_gn = gene_name
          best_known_pn = prot_name
        
        if species_id in REF_SPECIES:
          ref_known_gn = gene_name
          ref_known_pn = prot_name
          ref_id = pid
            
          if organism_db:
            key_proteins[organism_db].append(organism_db_xref)
            key_protein_ids.add(pid)
          
        cspecies.add(species)
      
      if (not best_known) or (not best_known_pn):
        if not best_known:
          best_known = ref_id or pid
          
        best_known_pn = ref_pn or prot_name
        best_known_gn = ref_gn or gene_name
      
      if key_proteins:
        key_protein_xrefs = []
        
        for organism_db in sorted(key_proteins):
          xrefs = ','.join(sorted(key_proteins[organism_db]))
          key_protein_xrefs.append(organism_db + ':' + xrefs)
        
        key_protein_xrefs = ';'.join(key_protein_xrefs)

      else:
        key_protein_xrefs = '-'
      
      key_protein_ids = ';'.join(key_protein_ids) or '-'
      
      line = f'{cid}\t{db_id}\t{name}\t{knownness:.1f}\t{num_prots}\t{len(cspecies)}\t{best_known}\t{best_known_gn}\t{best_known_pn}\t{key_protein_ids}\t{key_protein_xrefs}\n'
    
      write(line)
  
 
  cursor.close()
  connection.close()
  
  return file_path


def make_static_files(download_dir, static_dir, unknome_db_file_path=UNKNOME_DB_FILE_PATH):
  
  make_charts(download_dir, unknome_db_file_path) # For only current
  
  db_paths = glob.glob('unknome_db_*_*_*.sqlite') # All versions
  
  for db_path in db_paths:
    tsv_path = os.path.join(download_dir, TSV_CLUST_TABLE_FILE.format(_db_date_str(db_path)))
  
    if not os.path.exists(tsv_path):
      make_tsv_files(download_dir, db_path)


def make_tsv_files(download_dir, unknome_db_file_path=UNKNOME_DB_FILE_PATH):

  file_path = _make_summary_protein_table(download_dir, unknome_db_file_path)
  subprocess.call(['gzip', '-f', '-k', file_path])  
  
  file_path = _make_summary_cluster_table(download_dir, unknome_db_file_path)
  subprocess.call(['gzip', '-f', '-k', file_path])  
  

def make_charts(static_dir, unknome_db_file_path=UNKNOME_DB_FILE_PATH):

  save_path = os.path.join(static_dir, 'plot_knownness_orthologues.png')
  plot_knownness_vs_orthologues(unknome_db_file_path, save_path=save_path)
 
  save_path = os.path.join(static_dir, 'plot_knownness_change_stack.png')
  plot_knownness_change_stack(unknome_db_file_path, save_path=save_path)
 
  save_path = os.path.join(static_dir, 'plot_knownness_gains.png')
  plot_knownness_gains(unknome_db_file_path, plot_go_count=True, save_path=save_path)

  save_path = os.path.join(static_dir, 'plot_knownness_distribution.png')
  plot_knownness_distribution(unknome_db_file_path, save_path=save_path)


if __name__ == '__main__':
  
  download_dir = os.path.join(BASE_DIR, 'downloads')

  static_dir = os.path.join(BASE_DIR, 'static')
  
  make_static_files(download_dir, static_dir, UNKNOME_DB_FILE_PATH)
  
