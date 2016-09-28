# RBP_Motif_Search

A tool for browsing a database of RNA binding proteins (RBPs)

### RNA Binding Proteins

Genes are sequences encoded in DNA. These sequences are transcribed to short term RNA copies. The RNA molecules can encode proteins, regulate other genes or perform catalytic functions in the cell. Compared to DNA, they are short lived and are ultimately degraded.

Throughout the lifetime of an RNA it is always accompanied by proteins. Often these proteins recognize a specific sequence of nucleotides to bind to, called a motif

### RBP Motifs

The motifs that proteins bind to are often degenerate, meaning that certain mismatches are tolerated. For example, an RBP may prefer

GGCUAA

But may also bind to 

GGCUUA
GGCUAU
CGCUAA

For this reason, motifs are often represented a position weight matrices.

### Experimentally determined motifs

There are various computational and experimental approaches for determining the PWMs of particular RBPs. One study determined many such PWMs in high throughput and published their results online.

[A compendium of RNA-binding motifs for decoding gene regulation](https://www.ncbi.nlm.nih.gov/pubmed/23846655)

Full Citation:

Ray D, Kazan H, Cook KB, Weirauch MT, Najafabadi HS, Li X, Gueroussov S, Albu M, Zheng H, Yang A, et al. A compendium of RNA-binding motifs for decoding gene regulation. Nature. 2013;499(7457):172–177. doi: 10.1038/nature12311.


### A motif database

This project contains a database of motif information from the above study along with a lightweight Python tool for scanning sequences of interest for motifs.

Usage:

python run_motif_search.py my_sequences.fa

The input file is in [fasta format](run_motif_search.py). The tool outputs a few different reports about which motifs seem enriched. It can potentially support additional customization, but the UI for that is under development.










