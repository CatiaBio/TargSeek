# BLAST Databases

Reference proteomes used to check whether predicted epitopes are also found in
other organisms (cross-reactivity, e.g. human or cow proteins).

## Download

The databases are downloaded automatically in stage 1 of the download pipeline.
To download only the databases:

```bash
snakemake -s Snakefile_download all_blast_databases --cores 4
```

BLAST (`blastp`, `makeblastdb`) is included in the conda environment.

## Choose databases

In `config/config_download.yaml`, set each database to `true` or `false`:

```yaml
blast_databases:
  human: true
  cow: true
  mouse: true
  salmon: false
  ...
```

| Group | Databases |
|---|---|
| Standard | `human`, `viral_proteins`, `fungal_proteins` |
| Livestock | `cow`, `goat`, `sheep`, `buffalo`, `horse`, `pig`, `chicken` |
| Companion animals | `cat` (`dog` is not on the UniProt FTP at the moment) |
| Model organisms | `mouse`, `rat`, `rabbit`, `zebrafish` |
| Food / plants | `salmon`, `arabidopsis`, `rice`, `soybean` |
| Bacteria | `ecoli`, `salmonella`, `staph` |

Most proteomes are 15–80 MB. All of them together need about 500 MB.

## Refresh

Databases older than `blast_database_refresh.refresh_interval_months` (default: 6)
are downloaded again. Set `force_refresh_all: true` to refresh all of them now.

## Output

```
data/blast_databases/
├── human_proteome.fasta (+ .phr .pin .psq ...)   # one set per database
├── download.log
└── .blast_databases_complete.sentinel            # marks the download as done
```

## Manual use

```bash
blastp -query epitopes.fasta -db data/blast_databases/human_proteome.fasta \
       -evalue 1e-3 -outfmt 6 -out human_hits.tsv
```

## Problems

- **Download fails:** check your internet connection and that the UniProt/NCBI FTP servers are reachable.
- **`makeblastdb` not found:** activate the environment (`conda activate targseek`).
- **Database ignored:** names in the config are case-sensitive and must match the table above.
