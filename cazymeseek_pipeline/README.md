# CAZymeSeek/dbCAN research workflow

## Scope and evidence

This is a two-layer CAZyme workflow, not a new biological prediction model.

1. **Annotation/quantification layer.** `run_dbcan` implements the procedure in *Carbohydrate-active enzyme annotation in microbiomes using dbCAN* (bioRxiv 2024.01.10.575125v1): Trim Galore/Kraken2, MEGAHIT, Prokka, three-method CAZyme annotation, CGC finding, dbCAN-sub substrate mapping, dbCAN-PUL homology, majority voting, BWA mapping, Bedtools counts and `dbcan_utils` TPM aggregation.
2. **Knowledge layer.** The supplied CAZymeSeek asset is four ChromaDB stores: CAZypedia concepts, CAZy family records, dbCAN-seq CGCs and structured substrate reactions. Its local `all-MiniLM-L6-v2` model creates 384-dimensional cosine retrieval embeddings. It retrieves evidence; it does not annotate sequences or validate substrates.

The EB05 project PDF supplies CAZypedia, CAZy/UniProt, CGC, dbCAN-sub and RAG project context. The dbCAN PDF supplies every executable wet/dry analysis rule. The A/B/C labels below are an explicit reporting convention because neither PDF defines such labels.

## Layout

```text
cazymeseek_pipeline/
├── src/cazymeseek_pipeline/    Python modules
├── config/                     environment, example YAML and sample sheet
├── demo/                       offline schema/plot smoke test
├── output/                     generated results, logs and figures
├── db/                         mount point for run_dbcan data
├── vector_db/                  mount point for supplied four-store asset
└── README.md
```

## Deployment

### Conda

```bash
cd /data/cazymeseek_pipeline
conda env create -f config/environment.yml
conda activate cazymeseek
export PYTHONPATH="$PWD/src"
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
```

The versions in `environment.yml` reproduce the dbCAN PDF Materials list: Python 3.9, MEGAHIT 1.2.9, Trim Galore 0.6.0, Kraken2 2.1.1, Prokka 1.4, BWA 0.7.17, Samtools 1.7, HMMER 3.3, DIAMOND 2.1.8, BLAST 2.14, Bedtools 2.27.1 and Seqkit 2.5.1. `run_dbcan` is 4.0.0.

### Linux packages plus Python

APT repositories do not reliably provide the paper's exact versions. Use APT only for system libraries, then install the version-pinned bioinformatics environment with Conda.

```bash
sudo apt-get update
sudo apt-get install -y build-essential wget git default-jre
conda env create -f config/environment.yml
conda activate cazymeseek
```

## Databases

Build the dbCAN data with the exact installed release and record its release/date in the run log. The dbCAN PDF S4 supplies the database setup procedure; the command interface may change between releases.

```bash
mkdir -p db/run_dbcan
run_dbcan database --db_dir "$PWD/db/run_dbcan"
kraken2-build --standard --db db/kraken2
```

Mount the existing project asset rather than re-embedding it:

```bash
ln -s /data/vector_databases vector_db/current
# Set paths.vector_db_dir to /data/cazymeseek_pipeline/vector_db/current
```

It must contain all four stores and `models/all-MiniLM-L6-v2`. The PDF does not provide a downloadable CAZymeSeek archive, so this repository intentionally does not fabricate a URL.

## Run

```bash
# One individual-assembly sample
PYTHONPATH=src python -m cazymeseek_pipeline.pipeline --config config/config.yaml --sample Wet2014
# Sequential batch; one sample uses all configured threads
PYTHONPATH=src python -m cazymeseek_pipeline.pipeline --config config/config.yaml
# Offline structural smoke test
PYTHONPATH=src python demo/test_demo.py
```

### Assembly routes

| Route | Documented use and limitation | Resource implication |
|---|---|---|
| Individual | Default implementation; practical for most projects and preserves per-sample contig context. | Run one 40-thread sample at a time on the PDF reference host. |
| Co-assembly | Pools reads; can recover low-abundance genes and more complete genes, but can confuse strain-heterogeneous de Bruijn graphs. | Higher runtime/RAM; protocol cautions against large pools (for example, >5 samples). |
| Assembly-free | Direct reads-to-CAZyDB family profiling; cannot make CGCs and has higher false-positive/lower precision in the cited comparison. | Faster; use only as a separately labelled family-level comparison. |

The executable wrapper rejects non-`individual` modes rather than silently changing their biological meaning.

## Outputs and evidence labels

Standardized records retain: `sample_id`, `contig_id`, `gene_id`, `protein_sequence_id`, `protein_cluster_id`, `protein_primary_family`, `domain_index`, `domain_start`, `domain_end`, `cazy_class`, `family_id`, `subfamily_id`, `EC`, `annotation_evidence`, `selected_method`, `annotation_conflict`, `CGC_id`, `CGC_gene_composition`, `CAZyme_substrate`, `CGC_substrate_PUL`, `CGC_substrate_vote`, `TPM`, `RPM`, `RPKM`, `confidence_tier`, `source_database`, and `version`.

- **A:** dbCAN-sub subfamily with EC, or CGC substrate supported by dbCAN-PUL homology. The PDF recommends dbCAN-PUL homology over majority voting.
- **B:** a documented family call with annotation evidence. The PDF preference when methods disagree is dbCAN HMM > dbCAN-sub/eCAMI > DIAMOND.
- **C:** majority-vote or knowledge retrieval only. It is a candidate hypothesis, never a validated substrate result.

No numeric biological confidence threshold is assigned: neither source PDF defines one. MiniLM semantic score is only a cosine retrieval score.

## Five corrections in this release

1. **Conflict priority:** `annotation.py` records all dbCAN/dbCAN-sub/DIAMOND family calls and resolves only `protein_primary_family` as `dbCAN HMM > dbCAN-sub/eCAMI > DIAMOND`, the order in dbCAN Box 6. `selected_method` and `annotation_conflict` make every override auditable; individual `family_id` values remain domain-specific.
2. **Multi-domain preservation:** HMM and dbCAN-sub domain rows are expanded, with `domain_index`, `domain_start`, `domain_end`, and `subfamily_id`; one protein can therefore generate multiple output rows. No “first domain only” reduction occurs.
3. **EC matching:** EC is copied only from the curated dbCAN-sub/`overview.txt` release output. An uncurated family/subfamily remains blank: the PDF does not support an invented universal GH/GT-to-EC dictionary.
4. **PUL routes unchanged:** `run_dbcan --cgc_substrate` remains intact. `substrate.out` retains `CGC_substrate_PUL` and `CGC_substrate_vote` separately; it neither re-votes candidates nor replaces the document's preference of dbCAN-PUL homology over majority voting.
5. **Optional redundancy-aware abundance:** `deduplication.enabled: true` applies the PDF Box 8 MMseqs2 example (`>0.95` identity and `>0.95` coverage) before CDS mapping and reports a `protein_cluster_id`. Default `false` preserves the official P13 all-CDS `dbcan_utils` aggregate pathway. Cluster TPM is explicitly a project reporting mode, not a claim that the PDF requires it for all abundance analyses.

```bash
# Validate the fixes without databases or Linux bioinformatics executables.
PYTHONPATH=src python demo/test_demo.py
```

## Visualizations

```bash
PYTHONPATH=src python -m cazymeseek_pipeline.visualize heatmap --table output/Wet2014/abundance/fam_substrate_abund.out --output output/figures/substrates
PYTHONPATH=src python -m cazymeseek_pipeline.visualize bar --table output/Wet2014/abundance/fam_abund.out --label Family --abundance TPM --output output/figures/families
PYTHONPATH=src python -m cazymeseek_pipeline.visualize cgc --table output/Wet2014/dbcan/cgc_standard.out --cgc-id 'k141_41392|CGC3' --output output/figures/cgc
```

The heatmap summarizes dbCAN-sub substrate TPM; bars summarize family/subfamily/EC TPM; CGC structure uses `cgc_standard.out`; synteny requires `PUL.out`, `cgc_standard.out`, the best dbCAN-PUL GFF, and optional contig read depth, as specified in P16. PNG is 300 dpi and PDF is vector output.

## Troubleshooting and 40 CPU/128 GB operation

| Issue | Protocol-grounded check |
|---|---|
| Decontamination fails | Inspect Kraken2 report first; identify a source, then align to its reference. Carter2023 needed no removal. |
| Fragmented assembly | Check MEGAHIT output/N50; do not assume MAG binning is required; co-assembly is a separate comparison due to strain heterogeneity. |
| `run_dbcan` no result | Remove a pre-existing output directory and rerun; verify FAA and GFF `ID` values match for CGC prediction. |
| No retrieval hit | Check four collection directories, model path, and `family_id`/`cazyme_families_base`; semantic retrieval has no source-defined cutoff. |
| Few CGCs | Confirm prokaryotic input, matched GFF/FAA identifiers and contig continuity; CGC logic is not validated for eukaryotes. |

The PDF's reference execution takes about 33 h on 40 CPUs/128 GB. Gene prediction is its longest stage; substrate prediction is next; assembly has highest expected memory use. Allocate the whole host to one individual sample, retain stage logs/checkpoints, and schedule samples sequentially. If throughput is required, reduce per-job threads and reserve memory before parallelizing; do not oversubscribe MEGAHIT.
