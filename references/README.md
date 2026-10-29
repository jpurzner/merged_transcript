# References

## Primary publication

This repository holds the merged granule-neuron-precursor (GNP) developmental
timecourse underlying:

> **Ezh2 Delays Activation of Differentiation Genes During Normal Cerebellar
> Granule Neuron Development and in Medulloblastoma.**
> James Purzner, Alexander S Brown, Teresa Purzner, Lauren Ellis, Sara Broski,
> Ulrike Litzenburger, Kaytlin Andrews, Aryaman Sharma, Xin Wang,
> Michael D Taylor, Yoon-Jae Cho, Margaret T Fuller, Matthew P Scott.
> *eLife* (2025), reviewed preprint. DOI: [10.7554/eLife.105339.1](https://doi.org/10.7554/eLife.105339.1).
> Preprint: bioRxiv [10.1101/2024.11.21.624171](https://doi.org/10.1101/2024.11.21.624171).

Read it online:
- eLife reviewed preprint: https://elifesciences.org/reviewed-preprints/105339
- bioRxiv (full text + PDF): https://www.biorxiv.org/content/10.1101/2024.11.21.624171v1

The full PDF is intentionally **not** committed to keep the repo lean — follow
the links above.

### Abstract (summary)

Ezh2 (the H3K27me3 methyltransferase of PRC2) controls the timing of granule
neuron differentiation — the lineage of origin for SHH-subtype medulloblastoma.
H3K27me3 marks more than half of the genes whose expression rises as granule
neurons differentiate. Conditional Ezh2 knockout in GNPs triggers early
differentiation, but only *after* cell-cycle exit. In medulloblastoma cells the
Ezh2 inhibitor UNC1999 induces neuronal differentiation only when combined with
the CDK4/6 inhibitor palbociclib, which forces cell-cycle exit.

Organism: mouse (with human SHH-MB validation). Tissue: postnatal/embryonic
cerebellum; GNPs and granule neurons.

## How this repo relates to the paper

The merged timecourse here (Scott RNA-seq E15–P14 + Hatten TRAP array P0–P56)
is the data behind the paper's developmental gene-program classification. The
paper defines three gene sets by time of peak expression, reproduced exactly by
`R/figures/sankey_classification.R` (Figure S1):

| Program (paper) | n (paper) | This repo |
|-----------------|-----------|-----------|
| embryonic genes       | 783 | `gene_program == "embryonic"`      (peak E15), `fc_class == ">2 log2FC"` |
| proliferation genes   | 440 | `gene_program == "proliferation"`  (peak P0–P7, e.g. *Gli1*, *Ccnd1*), `fc_class == ">2 log2FC"` |
| differentiation genes | 787 | `gene_program == "differentiation"`(peak P14–P56, e.g. *Cacna1a*, *Cbx7*), `fc_class == ">2 log2FC"` |

The `gene_program` column in `data/timecourse_query_index.tsv` carries the
program by peak timing for every gene; the strict published sets additionally
require `fc_class == ">2 log2FC"` plus the scRNA-seq lineage / cell-cycle
exclusions applied in the Sankey (see `R/figures/sankey_classification.R` and
`analysis/3_gene_classification.Rmd`).

## Licensing

The paper itself is linked, not redistributed here. Reuse of the manuscript is
governed by the license on bioRxiv / the CC-BY terms of the eLife reviewed
preprint; cite the work above.

## BibTeX

```bibtex
@article{Purzner2025Ezh2,
  title   = {Ezh2 Delays Activation of Differentiation Genes During Normal
             Cerebellar Granule Neuron Development and in Medulloblastoma},
  author  = {Purzner, James and Brown, Alexander S and Purzner, Teresa and
             Ellis, Lauren and Broski, Sara and Litzenburger, Ulrike and
             Andrews, Kaytlin and Sharma, Aryaman and Wang, Xin and
             Taylor, Michael D and Cho, Yoon-Jae and Fuller, Margaret T and
             Scott, Matthew P},
  journal = {eLife},
  year    = {2025},
  note    = {Reviewed preprint},
  doi     = {10.7554/eLife.105339.1}
}
```
