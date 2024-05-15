---
aside: false
---

# Effects of mutations on pseudovirus entry in 293T cells 

## Interactive plot of mutation effects
The plot below shows how mutations affect cell entry, and allows you to zoom and mouseover sites and mutations. Positive values indicate better entry than the parental GPC while negative values indicate impaired entry. 

Here is an explanation of the key plot elements:
 - The zoom bar at the top of the plot shows different regions of GPC, and can be used to zoom in on specific sites.
 - The line plot summarizes the effects of mutations at each site, using the site smmary statistic specified by the interactive option at the bottom of the plot (eg, sum or mean effect of mutations at a site).
  - The heatmap shows how each individual mutation affects cell entry. The `x`'s indicate the amino-acid identity in the parental GPC strain, light gray indicates mutations that were not measured, and dark gray indicates mutations that are filtered out by one of the interactive sliders at the bottom of the plot. If you mouseover points, you can also see the measurement in each individual library, which provides a sense of the reproducibility between experiments. 
  - The options at the bottom of the plot allow you to adjust what is displayed. For instance, the *times_seen* slider filters for only mutations seen in at least that many unique variants per-library; larger values will generally give more accurate measurements. Similarly, you can filter by how many different experimental selections the mutation was measured in. You can also select to only show mutations accessible by some number of mutations to the codon in the parental sequence (eg, only accessible by single-nucleotide mutations).

Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/293T_entry_func_effects.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on entry into 293T cells">
    <Altair :showShadow="true" :spec-url="'htmls/293T_entry_func_effects.html'"></Altair>
</Figure>

## Numerical values of mutation effects on cell entry
The pre-filtered numerical data plotted on this page can be found [here](https://github.com/dms-vep/LASV_Josiah_GP_DMS/blob/main/results/filtered_func_effect_CSVs/293T_filtered_func_effects.csv). Alternatively, unfiltered data are available in [this CSV file](https://github.com/dms-vep/LASV_Josiah_GP_DMS/blob/main/results/func_effects/averages/293T_entry_func_effects.csv).
Note that these values have **not** been filtered by QC criteria like the *times seen*, so either make sure you understand the filters in the file or use the pre-filtered data.