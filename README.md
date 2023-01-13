## NextFlow Pipeline for discovering SARS-CoV-2 prolonged infection candidates in a variety of data sources

More documentation to come!

For now, here's an rather exhaustive, early example of how to run the workflow:
```
nextflow run long-infection-finder.nf \
--gisaid_metadata_dir /path/to/gisaid/metadata/directory \
--search_gisaid_metadata true \
--search_genbank_metadata true \
--search_genbank_seqs false \
--search_gisaid_seqs false \
--search_local_metadata false \
--min_date 2022-01-01 \
--max_date 2023-01-01 \
--geography Wisconsin \
--days_of_infection 240 \
--duration_of_interest 600
```

This command takes advantage of the many parameters (workflow settings) in the `nextflow.config` configuration file. These settings, in the order invoked above, include:
- Find the downloaded and unzipped GISAID metadata file, in TSV format, at the filepath specified after the flag `--gisaid_metadata_dir`. This path can be absolute or relative, and should not include the filename itself.
- Search the above-specified GISAID metadata for sequences classified as very old lineages by setting the `--search_gisaid_metadata` flag to `true`.
- Download and search metadata from NCBI GenBank for sequences with very old lineages by setting the `--search_genbank_metadata` flag to `true`.
- Do not pull and reclassify GenBank sequences in search of long infections, as set with `--search_genbank_seqs false`.
- Do not reclassify and search GISAID sequences for long infections, as set with `--search_gisaid_seqs false`.
- Do not search metadata for a locally available database of sequences, as set with `--search_local_metadata false`.
- Only search sequences in 2022 by setting `--min_date` to 2022-01-01 and `--max_date` to 2023-01-01. Dates *must* be formatted in this way.
- Only search for sequences collected in Wisconsin by setting `--geography` to "Wisconsin".
- Set a cutoff of ~8 months, 240 days, for identifying "anachronistic sequences," which are indicative of prolonged SARS-CoV-2 infection. These are sequecnes that classify as Pango lineages that arose ≥ 8 months ago, as specified with `--days_of_infection 240`.
- Only download sequences for infections that are 600+ days long, as set with `--duration_of_interest 600`.

The same command, pulling directly from GitHub, would be:
```
nextflow run nrminor/long-infection-finder.nf \
-latest \
--gisaid_metadata_dir /path/to/gisaid/metadata/directory \
--search_gisaid_metadata true \
--search_genbank_metadata true \
--search_genbank_seqs false \
--search_gisaid_seqs false \
--search_local_metadata false \
--min_date 2022-01-01 \
--max_date 2023-01-01 \
--geography Wisconsin \
--days_of_infection 240 \
--duration_of_interest 600
```