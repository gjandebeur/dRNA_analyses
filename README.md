# dRNA_analyses
Workflow for analyses of m6A Modifications Post Modkit Filtering


load dependencies first
```
library(data.table)
library(tidyverse) 
library(tidytable)
library(patchwork)
library(ggplot2)
library(GenomicRanges)
library(txdbmaker)
```

setwd and import m6A_subset from "m6A_subset_all.tsv"

### Plot data from filtered file, coverage & percent modified
This can and should be scaled up with more modification types and samples 
```
ggplot(m6A_subset, aes(x = avg_n_valid_cov, y = avg_percent_modified, color = mod)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, linewidth = 1.2) +
  scale_x_log10(limits = c(10, 10000)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Valid Coverage (log10)",
    y = "Percent Modified",
    color = "Modification",
    title = "Coverage and Percent Modification"
  ) +
  theme_classic(base_size = 14)
```
<img width="1200" height="800" alt="githubexample_coverage_percentmodified" src="https://github.com/user-attachments/assets/9aa987d0-7218-4c6d-8df4-c3c9030b8554" />


load in txdb from gtf file (provided in files)

```
txdb <- makeTxDbFromGFF("at_ensembl_plants.gtf")
k <- keys(txdb, keytype = "TXNAME")  
tx2gene_at <- AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
```

intersect txdb with the m6A_subset dataframe

```
#intersecting files
m6A_subset <- m6A_subset %>%
  inner_join(., tx2gene_at, by=c("tx"= "TXNAME")) %>%
  unite(position, chrom:end, remove =F) 

#changing format to make reproducible
m6A_subset <- m6A_subset %>%
  separate(position, into = c("chrom", "start", "end"), sep = "_", remove = FALSE) %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  )
```

### Now running Genomic Ranges to split tx into their genomic regions
```
mods <- GRanges(
  seqnames = m6A_subset$chrom,
  ranges = IRanges(start = m6A_subset$start, end = m6A_subset$end),
  mod = m6A_subset$mod,
  gene_id = m6A_subset$gene_id
)

cds_gr <- cds(txdb)
five_utr_gr <- fiveUTRsByTranscript(txdb, use.names = TRUE) |> unlist()
three_utr_gr <- threeUTRsByTranscript(txdb, use.names = TRUE) |> unlist()
introns_gr <- intronsByTranscript(txdb, use.names = TRUE) |> unlist()
```

### Now running function to find each region's hits
```
assign_region <- function(mods, cds, five_utr, three_utr, introns) {
  region <- rep("intergenic", length(mods))  # default assignment

  hits_cds <- findOverlaps(mods, cds)
  region[queryHits(hits_cds)] <- "CDS"
  hits_5utr <- findOverlaps(mods, five_utr)
  idx_5utr <- setdiff(queryHits(hits_5utr), queryHits(hits_cds))
  region[idx_5utr] <- "5'UTR"  
  hits_3utr <- findOverlaps(mods, three_utr)
  idx_3utr <- setdiff(queryHits(hits_3utr), c(queryHits(hits_cds), idx_5utr))
  region[idx_3utr] <- "3'UTR"
  
  hits_introns <- findOverlaps(mods, introns)
  idx_introns <- setdiff(queryHits(hits_introns), c(queryHits(hits_cds), idx_5utr, idx_3utr))
  region[idx_introns] <- "intron"
  
  return(region)
}

# add region tag to the original dataframe for downstream analyses
m6A_subset$region <- assign_region(mods_gr, cds_gr, five_utr_gr, three_utr_gr, introns_gr)

#Lastly, filter by descending groups 
summary_filtered_data_byregion <-  m6A_subset %>%
  group_by(region, mod) %>%
  summarise(n_mods = n(), .groups = "drop") %>%
  arrange(desc(n_mods))

# For overview of graph
table(m6A_subset$region)
```

### Now plotting genomic ranges
```
region_mod_summary <- m6A_subset %>%
  group_by(region, mod) %>%
  summarise(n_mods = n(), .groups = "drop")

filtered_region_mod_summary <- region_mod_summary %>%
  filter(mod %in% c("pseU", "m6A"))

#Now Actually Plot
## This can be easily scaled up with other modification types or species

ggplot(filtered_region_mod_summary, aes(x = region, y = n_mods, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Modification Counts per Region (pseU & m6A)",
    x = "Genomic Region",
    y = NULL,
    fill = "Genomic Region"
  ) +
  scale_fill_aaas() +
  facet_wrap(~ mod, scales = "free_y") +
  theme_minimal(24) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank()  # this removes tick labels under bars
  )
```
<img width="1200" height="900" alt="githubexample_genomicregions" src="https://github.com/user-attachments/assets/fd1854b8-5856-433d-ab47-8c1002879210" />


### Next export the subset file as bed to run bedtools intersect (must be done unix) to line up with reference file and output as fasta

```
m6A_subset %>%
  filter(code == "a") %>%
  dplyr::select(chrom, start, end, code, n_valid_cov, strand) %>%
  distinct(chrom, start, end, .keep_all = T) %>%
  arrange(chrom, start) %>%
  mutate(start = start - 3,
         end = end + 3) %>%
  write.table("C:/extract/to/bedfile.bed",
              col.names = F, row.names = F, sep = "\t", quote = F)
```

### Next code can NOT be run in Rstudio, must export to a UNIX shell

next part must be run on Unix with bedtools installed to run the intersect on the bed file with the known fasta reference
```
./bedtools getfasta -s -fi "reference_toplevel.fa" -bed "C:/extract/to/bedfile.bed" > "bedfile_nowfa.fasta"

samtools faidx bedfile_nowfa.fasta"
```

### Using pA lengths from original basecalling for analyses

```
m6A_subset_tails <- m6A_subset %>%
  mutate(
    tail_length_group = case_when(
      mean_tail_length <= quantile(mean_tail_length, 0.25, na.rm = TRUE) ~ "short",
      mean_tail_length >= quantile(mean_tail_length, 0.75, na.rm = TRUE) ~ "long",
      TRUE ~ "medium"
    )
  )
```

### Label PolyA tail Categories for Visualization 



