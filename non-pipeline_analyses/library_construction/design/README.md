# Library Design

Kate started with a download of all 1183 mammarenavirus GP *protein* sequences available on genbank on 11/19/20 as `mammarenavirus_gpcs_all.fasta`.

Exact search in the protein database:
    
    `(("Mammarenavirus"[Organism] OR mammarenavirus[All Fields]) AND gpc[All Fields]) AND viruses[filter]`

This did **not** filter out partial cds. 

Kate then did subsequent searches of GenBank and the VIPR database to look for more complete sets of sequences. 

These downloaoded sequences are in the `seq_downloads` directory. These searches are compared in the `design.ipynb` notebook.

The `Mammarenavirus glycoprotein` GenBank search (same search terms as above, but with `glycoprotein` spelled out instead of `gpc`) yielded the most full-length glycoprotein sequences from New World and Old World mammarenaviruses.    

The `design.ipynb` notebook will filters the sequence fasta files, creates multiple sequence alignments of all glycoprotein sequences and just LASV GPC sequences, constructs phylogenetic trees, and further analyzes these alignments for library design.


## Main Goal: 

Although we are recovering high titers of LASV GP-pseudotyped lentivirus from cells transduced with WT LASV GP, these titers drop more than expected when only 1/3 of the library is functional.

In order to rescue a library with the greatest diversity of functional variants as possible, I wanted to look into using natural sequences to identify regions of LASV GP to *exclude* from mutagenesis.
Regions with complete conservation across all mammarenavirus GPC sequences are potential candidates for exclusion from mutagenesis as mutagenizing these regions may be more likely to result in non-functional proteins. 


### Main Findings:

There are 260 sites fully conserved in the LASV GP alignment. Excluding those sites would more than halve the library, but I want to identify mutations that may not have been sampled naturally. 

If we include all Old World GPCs, there are 102 sites that are conserved. However, the alignment is mostly made up of LASV sequences (706/1050), so I worry I am not really fully sampling the diversity of Old World GPCs (as many of those 344 non-LASV sequences are also New World viruses). 
As such, I am hesitant to limit the library only to sites with natural diversity. 

Given recent high titers of recovered LASV GP Josiah (re-optimized) pseudotyped lentiviruses (>1e5 TU/mL), I am less concerned about large portions of the library being non-functional.
Therefore, I am **not** going to exclude conserved sites from the library. 

I am also **not** going to exclude the signal peptide, since it is actually included as a final component of the LASV GPC and may have functional roles. However, I am going to exclude mutating the transmembrane domain and cytoplasmic tail. Mutations in the transmembrane domain are likely going to be detrimental and, given, the fairly large amount of natural diversity in the cytoplasmic tail, mutations in the tail are quite likely to either be fairly artificial artifacts of the pseudotyped system or simply neutral.


## Side Goal:

* Examine and get a better handle on mammarenavirus GP diversity. 

    These analyses clearly show different amounts of divergence from LASV GP Josiah for New World and Old World mammarenavirus GPCs.

    Additional phylogenetic analyses (including better annotation of branch length) would be necessary to fully confirm this, but mammarenavirus GPC evolution does  seem to be more geographically, rather than temporally clustered (i.e. the tree is not ladder-like). This is consistent with spillover from geographically isolated rodent populations.


* Perhaps add in additional sequences from other arenavirus species to better understand arenavrius GP diversity.

    Other arenavirus sequences are from non-mammalian hosts and I did not think it was relevant to add those.

## Conclusions/Next Steps

Given the analysis of currently available Mammarenavirus GPC sequences contained in the `design.ipynb` notebook, I decided to not limit mutation sites for constructing a LASV GP library.
I was initially going to not mutate the TM or cytoplasmic tail, but found a couple of papers ([Zhang, et al., 2019 ](https://link.springer.com/article/10.1007/s13238-018-0604-x?utm_source=TrendMD&utm_medium=cpc&utm_campaign=Protein_%2526_Cell_TrendMD_1&origin=03d9293ceeea6f4df236242eb107bb3e) and [Willard, et al., 2019](https://www.mdpi.com/2076-0817/8/1/1/htm)) that look at mutating the TM and cytoplasmic tail, so I decided to mutate the entire GP. 

It may also be worth looking at the regions of conserved sites on the structure and seeing if it generally correlates with RSA, but I'm not sure how much it is worth my time to do this versus get a library ordered. 
However, I briefly did look to see if these sites were typically located in the core or surface of GP. 
There was little pattern for LASV GP conserved residues, the OW and all mammarenavirus conserved residues were potentially more commonly found in the core/GP2, but there were not obvious patches to avoid mutating.

A typical codon tiling library spanning all regions of LASV GP using the LASV GP Josiah (re-codon-optimized) sequence is designed in the `codon_tiling_primers` subdirectory of the `library_design` directory.


## Ordering primers

IDT's oPools (oligo pools) only allow for `N` or `K` ambiguous nucleotides, so the`S` nucleotide cannot be used.
In order to have the equivalent of `NNS`, `NNG` and `NNC` primers were designed (`CNN` and `GNN` in reverse complement). 

These primers were designed by using a slightly modified script in the ([CodonTilingPrimers](https://github.com/jbloomlab/CodonTilingPrimers)) repo to allow for ambiguous codons beyond `NNN`, `NNS`, or `NNK`. 
By allowing for `NNG` and `NNC` codons to be designed, the script can simply be run twice (once for each codon) and the `NNG`/`NNC` forward primers combined for one pool and the `CNN`/`GNN` primers combined for the other pool.

The exact script for making primers was copied to this repo in the `codon_tiling_primers` directory as `2021Jan_create_primers.py`

Based on Tm and primer length analyses in the different Tm excel files (graphs in Tm61-62 file), I decided to change the Tm range for primer design to 60.5-61.5.

To make libraries, the following commands were run from within the `codon_tiling_primer` directory. 

```
python 2021Jan_create_primers.py 210112_lasvgpc_josiah_reopt_excessflank.txt LASV_GP_Josiah_reopt 1 210113_LASVGP_Josiah_reopt_Primers_NNG.csv --ambiguous_codon NNG --output opools --minprimertm 60.5 --maxprimertm 61.5
```
and
```
python 2021Jan_create_primers.py 210112_lasvgpc_josiah_reopt_excessflank.txt LASV_GP_Josiah_reopt 1 210113_LASVGP_Josiah_reopt_Primers_NNC.csv --ambiguous_codon NNC --output opools --minprimertm 60.5 --maxprimertm 61.5
```

The resulting forward and reverse `NNG` and `NNC` (rev: `CNN`, `GNN`) primer files were combined (via copy/pasting) into one excel spreadsheet for ordering opools.

The ordered primers are in the `20210114_LASVGP_Josiah_reopt_AllPrimers.xlsx` file in the `codon_tiling_primers` directory.
