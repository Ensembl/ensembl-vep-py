---
myst:
  html_meta:
    "description lang=en": |
      The input and output file formats Ensembl VEP reads and writes.
---

# Data formats

The files Ensembl VEP reads and writes. See [Choosing annotation sources](../user_guide/how_to/annotation_sources.md) and [Writing output and run summaries](../user_guide/how_to/output.md) for how to use them.

## Input

### Variants

Variants are read from a [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) file, given with `--input_file`. The file can be plain text or compressed with `bgzip`. A file without header lines is accepted.

Variants should be sorted by position within each sequence region. A run refuses input that is clearly unsorted, rather than annotating it with transcripts from the wrong region.

Symbolic alleles such as `<DEL>`, breakend notation, and records with an `SVTYPE` INFO field are read as structural variants; all other records are read as short variants.

### Transcripts

Transcripts are read from one [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) file (`--gff`) or one [GTF](https://www.ensembl.org/info/website/upload/gff.html) file (`--gtf`). The file can be plain text or compressed with `bgzip`, and can have a `tabix` index.

A transcript is assembled from its own record, its exon and CDS records, and its parent gene record. With a `tabix` index, these are read for the regions around the variants. Without an index, the file is read from the start, and all the records of a gene must be next to each other.

### Reference sequence

The reference sequence is read from a [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) file given with `--fasta`. The file can be plain text or compressed with `bgzip`. An index is optional: a plain file uses a `.fai` index, and a `bgzip` compressed file uses both a `.fai` and a `.gzi` index.

### Sequence-region synonyms

A synonyms file, given with `--synonyms`, is plain text with two whitespace-separated columns: the name used by the reference sequence and a synonym for it. Lines with fewer than two names are skipped; columns after the second are ignored with a warning.

```text
22	NC_000022.11
X	NC_000023.11
```

### Configuration files

`--config` reads a configuration file in TOML or in the Perl VEP format. See [Configuring a run](../user_guide/how_to/configuring_a_run.md#configuration-files).

## Output

### VCF output

The default output format, `--output_format vcf`, writes each input record with its annotations added in a `CSQ` INFO field. The output header adds three lines to the input header:

- `##VEP`, with the version and the time of the run;
- `##INFO=<ID=CSQ,...>`, which declares the fields of the `CSQ` value;
- `##VEP-command-line`, with the command that produced the file.

The `CSQ` value holds one comma-separated entry per allele and feature. Each entry holds these fields, separated by `|`:

| Field | Content |
| --- | --- |
| `Allele` | The variant allele the entry describes. |
| `Consequence` | The consequence terms, separated by `&`. |
| `IMPACT` | The impact rating of the most severe term: `HIGH`, `MODERATE`, `LOW` or `MODIFIER`. |
| `SYMBOL` | The gene symbol. |
| `Gene` | The gene stable ID. |
| `Feature_type` | The type of the feature, for example `Transcript`. |
| `Feature` | The feature stable ID. |
| `BIOTYPE` | The biotype of the transcript. |
| `EXON` | The exon number, out of the total number of exons. |
| `INTRON` | The intron number, out of the total number of introns. |
| `HGVSc` | The HGVS coding sequence name. |
| `HGVSp` | The HGVS protein sequence name. |
| `cDNA_position` | The position in the cDNA. |
| `CDS_position` | The position in the coding sequence. |
| `Protein_position` | The position in the protein. |
| `Amino_acids` | The reference and variant amino acids. |
| `Codons` | The reference and variant codons, with the changed bases in upper case. |
| `Existing_variation` | Identifiers of known variants at the same position. |

The fields follow the layout Perl VEP writes, so tools that read its `CSQ` field can read this output. Fields the current implementation does not fill yet are left empty.

### JSON output

`--output_format json` writes one JSON object per line, one line per variant, in the Ensembl VEP JSON format. A variant object holds:

- `id`, `seq_region_name`, `start`, `end`, `strand` and `allele_string`, describing the variant;
- `input`, the input record it was read from;
- `most_severe_consequence` and `variant_class`;
- `transcript_consequences`, a list with one object per allele and transcript. Each object holds `variant_allele`, `consequence_terms`, `impact`, `transcript_id`, `gene_id` and `strand`, and, where they apply, the cDNA, CDS and protein positions, `amino_acids` and `codons`.

### Run summary

`vep run` writes a run summary beside its output, and `vep stats` builds one from an annotated VCF file. `--stats_format` (for `vep run`) or `--output_format` (for `vep stats`) selects the format:

- `html`, a page with tables and charts;
- `text`, plain-text tables;
- `json`, the same counts as a JSON document.
