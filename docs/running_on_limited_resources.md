# Running on limited resources

This machine (WSL2) has **12 cores but only ~7.7 GB of RAM**, which WSL takes as
about half the host's memory. RAM, not CPU, is the limit — a run was killed when a
PyTorch-based tool was started while the download pipeline was already going.

## Rules of thumb

1. **Run one heavy job at a time.** Never start ConSurf, BepiPred, DeepLocPro or
   DiscoTope while a download or analysis run is going.
2. **Cap the thread-hungry libraries.** NumPy, PyTorch and XGBoost grab all 12 cores
   by default, each thread costing memory.
3. **Use `--cores 2`** for the pipelines. Most download steps wait on APIs, so extra
   cores buy nothing and cost memory.
4. **Be nice to the desktop**: `nice -n 10` keeps VS Code responsive.

## The command to use

```bash
cd ~/Projects/TargSeek
OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 \
  nice -n 10 conda run --no-capture-output -n targseek \
  snakemake -s Snakefile_download all_download_data --cores 2 --keep-going \
  > logs/download_$(date +%Y%m%d_%H%M).log 2>&1
```

Same pattern for `Snakefile_analysis all_analysis`.

`--keep-going` matters: one failing step (a flaky API, say) then does not stop the
rest of the run.

## Watching it

```bash
# memory
free -m | awk 'NR==2{printf "used %dMB / %dMB (%.0f%%)\n",$3,$2,$3/$2*100}'

# biggest processes
ps -eo pmem,rss,comm --sort=-rss | head -5

# progress
tail -f logs/<the log file>
```

If "used" gets near 7000 MB, stop the run rather than letting the kernel kill it.

## Memory cost per step

| Step | Cost | Note |
|---|---|---|
| BacDive classification | low | API-bound; the cache file is ~150 MB |
| QuickGO fetch | low–medium | grows with the annotation set |
| Gene aliases | low | API-bound, resumable via cache |
| Coverage assessment | medium | loads the gene–species cache |
| Sequence / structure download | low | network-bound |
| MAFFT | medium | scales with sequences per gene |
| ConSurf (rate4site) | medium | one gene at a time |
| BepiPred 3.0 | **high** | ESM-2 model, run alone |
| DeepLocPro | **high** | ESM-2 model, run alone |
| DiscoTope-3.0 | **high** | ESM-IF1 model, run alone, `cpu_only: true` |

For the three high-cost tools, also set the thread caps above and give them the
machine to themselves.

## Raising the WSL memory limit (optional)

WSL currently gets ~7.7 GB. If the host has 16 GB you can give WSL more by creating
`C:\Users\catia\.wslconfig` **on the Windows side**:

```ini
[wsl2]
memory=10GB
swap=4GB
```

Then run `wsl --shutdown` in PowerShell and reopen the terminal. Leave the host at
least 4–6 GB for Windows itself.

## Interrupted runs

Nothing is lost. The caches under `cache/` and the pipeline's sentinel files mean a
re-run skips finished work and resumes partial downloads. Just repeat the same
command.
