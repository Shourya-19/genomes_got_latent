# SLIME (Similarity-based Latent gene Identification MEthod)
Pipeline to analyse any genome and predict latent sequences
## Download and setup the repository
```
git clone your_repo
conda env create -f environment.yml
conda activate slime
```
slime run config.yaml

## 1. Translating the genome in all six reading frames
Translate the genome in all reading frames. The default length cutoff is 60 if not specified.
``` bash
python translation.py test_genome.fasta all_translations.txt
```
## 2. Make an organism-specific proteome database
Download the proteome in FASTA format. 
``` bash
makeblastdb -in organism_proteome.fasta -title organism_db -dbtype prot -out organism_db/organism_db 
```
 * The -in flag states the fasta file to create the database from
 * The -title flag gives the database a title
 * The -dbtype flag says whether it is protein (prot) or nucleotide (nucl)
 * The -out flag is for the name of the database
 * The -parse_seqids states we want to retain the full names of each sequence

## 3. Run organism-specific BLAST
The evalue is set to 1e-5, and `-num_threads` can be adjusted based on the system's configuration. **One does not need to provide the sequences individually.**
```
blastp -query all_translations.txt -db organism_db/organism_db -out organism_blast_output.tsv -evalue 1e-5 -num_threads 4
```
## 4. Extract only sequences that give no hits in step 3.
```
python filter_nohit.py all_translations.txt organism_blast_output.tsv no_hit_org_blast.fasta
