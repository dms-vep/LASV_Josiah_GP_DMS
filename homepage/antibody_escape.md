---
aside: true
---

# Escape from monoclonal antibodies

The plots below show how mutations affect neutralization by a panel of monoclonal antibodies (25.10C, 12.1F, 37.7H, 25.6A, 37.2D, and 8.9F). The plots are interactive, and allow you to zoom and mouseover sites and mutations. Positive values indicate reduced neutralization by antibody while negative values indicate more neutralization.

Here is an explanation of the key elements for each plot:
 - The zoom bar at the top of the plot shows different regions of GPC, and can be used to zoom in on specific sites.
 - The line plot summarizes the neutralization escape caused by mutations at each site (larger values indicate more escape). The escape at a site is quantified using the site summary statistic specified by the interactive option at the bottom of the plot (eg, sum or mean effect of mutations at a site).
 - The heatmaps then show how each individual mutation affects antibody neutralization. The `x`'s indicate the amino-acid identity in the parental GPC strain, light gray indicates mutations that were not measured, and dark gray indicates mutations that are filtered out by one of the interactive sliders at the by one of the interactive sliders at the bottom of the plot (eg, they have a very negative effect on cell entry). You can mouse over any of the points. You can mouse over any of the heatmap cells for details.
 - The options at the bottom of the plot allow you to interactively adjust what is displayed. For instance, 
   the *minimum mutation entry in 293T cells* only shows mutations with at least some minimal cell entry
   function (and grays out mutations that are more deleterious). You can also select to floor the escape at zero to show / not show "negative" escape values. 


## Antibody 25.10C
- The plot below shows how mutations affect neutralization by the monoclonal antibody 25.10C. Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/2510C_mut_effect.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on antibody 25.10C escape">
    <Altair :showShadow="true" :spec-url="'htmls/2510C_mut_effect.html'"></Altair>
</Figure>

## Antibody 12.1F
- The plot below shows how mutations affect neutralization by the monoclonal antibody 12.1F. Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/121F_mut_effect.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on antibody 12.1F escape">
    <Altair :showShadow="true" :spec-url="'htmls/121F_mut_effect.html'"></Altair>
</Figure>

## Antibody 37.7H
- The plot below shows how mutations affect neutralization by the monoclonal antibody 37.7H. Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/377H_mut_effect.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on antibody 37.7H escape">
    <Altair :showShadow="true" :spec-url="'htmls/377H_mut_effect.html'"></Altair>
</Figure>

## Antibody 25.6A
- The plot below shows how mutations affect neutralization by the monoclonal antibody 25.6A. Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/256A_mut_effect.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on antibody 25.6A escape">
    <Altair :showShadow="true" :spec-url="'htmls/256A_mut_effect.html'"></Altair>
</Figure>

## Antibody 37.2D
- The plot below shows how mutations affect neutralization by the monoclonal antibody 37.2D. Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/372D_mut_effect.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on antibody 37.2D escape">
    <Altair :showShadow="true" :spec-url="'htmls/372D_mut_effect.html'"></Altair>
</Figure>

## Antibody 8.9F
- The plot below shows how mutations affect neutralization by the monoclonal antibody 8.9F. Click on the expansion box in the upper right of the plot to enlarge it for easier viewing, or [click here](/htmls/89F_mut_effect.html){target="_self"} to open the plot in a stand-alone window.

<Figure caption="Interactive plot showing effects of mutations on antibody 8.9F escape">
    <Altair :showShadow="true" :spec-url="'htmls/89F_mut_effect.html'"></Altair>
</Figure>

## Numerical values of escape

For per-antibody escape values, see [these CSV files](https://github.com/dms-vep/LASV_Josiah_GP_DMS/tree/main/results/filtered_antibody_escape_CSVs).
