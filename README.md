# SLIME (Similarity-based Latent gene Identification MEthod)
Pipeline to analyse any genome and predict latent sequences

## 1. Translating the genome in all six reading frames
Translate the genome in all reading frames. The default length cutoff is 60 if not specified.
``` bash
python 6_Frame_Translation.py input.fasta all_translations.txt
