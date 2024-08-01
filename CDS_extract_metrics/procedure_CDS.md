# Pipeline for Codon Metrics Analysis

This pipeline describes the steps to process genomic annotations, extract coding sequences (CDS), filter them, and compute various codon usage and sequence-based measures.

## Genome Annotation and CDS Extraction
### Data Sources and Conversion

**RIKEN:**

Convert GFF3 and FASTA files to CDS FASTA files:

```bash
gffread file.gff3 -g file.fa -x cds.fa
```
**JGI:**

_For GFF3 files:_
```bash
gffread file.gff3 -g file.fa -x cds.fa
```
_For GFF2 files:_

```bash
python EXTRACT_CDS_GFF2.py input.gff2 output.cds.fa
```
**iGenolevures:**

Convert EMBL files to FASTA:

```bash
python genbank_to_fasta_update.py input.embl output.cds.fa
```
**NCBI:**

Handle alternative splicing and convert:

```bash
python genbank_to_fasta_update.py input.gbff output.cds.fa
```

### Verification

Ensure CDS names are unique and valid:

```python
# Check for duplicates and invalid characters
python checkCDS.py input.cds.fa
```
## CDS Filtration

Filter CDS based on specified criteria:

```python
python filter_cds.py input.cds.fa filtered.cds.fa
```
## Codon Counting

Count codons for both normal and wobble positions:

```python
python Codon_count_V3.py input.cds.fa codon_counts.txt
```
## Sequence-Based Measures

Calculate nucleotide counts and GC content:

```python
python calculate_gc_content.py input.cds.fa gc_content.txt
```
Extract CDS coordinates and calculate surrounding GC content:

```python
python extract_coordinates.py input.gff3
python calculate_surrounding_gc.py input.fa coordinates.txt surrounding_gc.txt
```

## Codon Usage and Optimal Codons

Calculate various codon usage metrics:

```bash
# Determine optimal codons using codonW
codonw input.cds.fa -seq -all_indices -output codonw_output.txt
```
Compute CAI, CBI, and Fop:

```python
python calculate_optimal_codons.py codonw_output.txt optimal_codons.txt
```

## Codon Context Measures

Perform context census and calculate related metrics:

```python
python contextCensus.py input.cds.fa context_counts.txt
python calculate_context_measures.py context_counts.txt context_measures.txt
```
## Protein Composition from CDS

Calculate physico-chemical properties of proteins:

```python
python calculate_protein_properties.py input.cds.fa protein_properties.txt
```
## Aggregation of Measures

Aggregate all computed measures into final tables:

```python
python aggregate_measures.py filtered.cds.fa codon_counts.txt optimal_codons.txt context_measures.txt protein_properties.txt final_measures.txt
```
## Output Files

- Filtered CDS Sequences: One FASTA file per species.
- Codon Counts: Two files per species (complete and wobble).
- CodonW Analyses: One summary file per species.
- Context Counts: Two files per species (complete and wobble).
- Measures: Three files per species (complete, wobble, aminoacid).
- Taxonomy File: One file (genome-species.csv).

This pipeline ensures a comprehensive analysis of codon usage patterns and related metrics, providing valuable insights for taxonomic classification studies.