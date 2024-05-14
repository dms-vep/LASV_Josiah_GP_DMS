---
layout: home

hero:
  name: "Pseudovirus deep mutational scanning of Lassa virus glycoprotein complex"
  tagline: "Experimental measurements of how mutations to the glycoprotein complex (GPC) from the lineage IV Lassa virus Josiah strain affect cell entry and antibody escape."
  image: Main_website_image.png
features:
  - title: Summary
    details: Summary of how mutations affect cell entry and antibody escape
    link: /summary
  - title: Cell entry
    details: Effects of mutations on pseudovirus entry in 293T cells
    link: /cell_entry
  - title: antibody escape
    details: Effects of mutations on neutralization by monoclonal antibodies
    link: /antibody_escape
  - title: Interactive phylogenetic tree
    details: Nextstrain tree of Lassa GPC sequences that can be colored by antibody escape
    link: https://nextstrain.org/groups/dms-phenotype/lassa/lassa-GPC
  - title: Experiments and biosafety
    details: Explanation of pseudovirus deep mutational scanning
    link: /experiments_and_biosafety
  - title: Computational pipeline
    details: Details on the computational pipeline
    link: /pipeline_information
---

## About this site
This website contains interactive plots and numerical results from [pseudovirus deep mutational scanning](https://doi.org/10.1016/j.cell.2023.02.001) experiments that measure the effects of mutations to the glycoprotein complex (GPC) of the lineage IV Lassa virus Josiah strain on cell entry and antibody escape.
This study was led by Caleb Carr and Kate Crawford in the Bloom lab.

The links in the boxes above take you to interactive plots or descriptions of different aspects of the study.
For a high-level overview, see the [summary](summary){target="_self"} of how mutations affect cell entry and antibody escape.
To delve into the data in more detail, click on the boxes above for each individual phenotype.

You can also examine the output of the full [computational pipeline](pipeline_information){target="_self"} and look at the underlying code [on GitHub](https://github.com/dms-vep/LASV_Josiah_GP_DMS.git).
Here are CSV files of the numerical values pre-filtered for high-confidence values for [cell entry](https://github.com/dms-vep/LASV_Josiah_GP_DMS/blob/main/results/filtered_func_effect_CSVs/293T_filtered_func_effects.csv) and [antibody escape](https://github.com/dms-vep/LASV_Josiah_GP_DMS/tree/main/results/filtered_antibody_escape_CSVs).

Note the experiments use [single-cycle replicative pseudoviruses that can be safely studied at biosafety-level 2](experiments_and_biosafety){target="_self"}.

