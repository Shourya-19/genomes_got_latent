# SLIME (Similarity-based Latent gene Identification MEthod)
Pipeline to analyse any genome and predict latent sequences
## Download and setup the repository
```
git clone https://github.com/Shourya-19/slime.git
conda env create -f environment.yml
conda activate slime
```
## 1. Translating the genome in all six reading frames
Translate the genome in all reading frames. The default length cutoff is 60 if not specified.
``` bash
python translation.py test_genome.fasta all_translations.txt
```
## 2. Make an proteome databases
Note: Download the proteome in FASTA format. 
### a. Organism specific database
``` bash
makeblastdb -in organism_proteome.fasta -title organism_db -dbtype prot -out organism_db/organism_db 
```
### b. UniRef90
The paper uses the UniRef90 dataset available here. https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/
Make a database using:
```
makeblastdb -in UniRef90.fasta -title UniRef90 -dbtype prot -out UniRef90_db/UniRef90_db 
```
## 3. Run organism-specific BLAST
The evalue is set to 1e-5, and `-num_threads` can be adjusted based on the system's configuration. **One does not need to provide the sequences individually.**
```
blastp -query all_translations.txt -db organism_db/organism_db -out organism_blast_output.tsv -evalue 1e-5 -num_threads 4 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
```
## 4. Extract only sequences that give no hits in step 3.
```
python filter_nohit.py all_translations.txt organism_blast_output.tsv no_hit_org_blast.fasta
```
## 5. Run UniRef90 BLAST on no_hits file obtained in step 4.
```
blastp -query no_hit_org_blast.fasta -db UniRef90_db/UniRef90_db  -out uniref_blast.tsv -evalue 1e-5 -num_threads 4 -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
```

## 6. Analyse and shortlist putative latent proteins
The user has to give a) the blast out file from step 5, b) translation file from step 1, c) output file name and d) their organisms Genus.
```
python analyse_uniref_hits.py uniref_blast.tsv all_translations.txt latent_proteins.txt genus
```
