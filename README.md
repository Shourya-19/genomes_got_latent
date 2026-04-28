# SLIME (Similarity-based Latent gene Identification MEthod)
Pipeline to analyse any genome and predict latent sequences

## 1. Translating the genome in all six reading frames
Translate the genome in all reading frames. The default length cutoff is 60 if not specified.
``` bash
python translation.py test_genome.fasta test_output_translations.txt
