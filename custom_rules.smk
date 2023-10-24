"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""

rule spatial_distances:
    """Get spatial distances from PDB."""
    input: 
        pdb="data/7puy.pdb",
    output:
        csv="results/spatial_distances/spatial_distances.csv",
    params:
        target_chains=["a", "b", "c", "A", "B", "C"],
    log:
        log="results/logs/spatial_distances.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    script:
        "scripts/spatial_distances.py"


rule averaged_ridgeplot_func_scores:
    """
    Average the functional score variants by type
    and display as a ridgeplot for the final functional 
    selections used in the averages.
    """
    input:
        func_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        nb="notebooks/visualize_func_scores_by_variant_type.ipynb",
    params:
        libA_1="LibA-220823-293T-1",
        libA_2="LibA-220823-293T-2",
        libA_3="LibA-220907-293T-1",
        libA_4="LibA-220907-293T-2",
        libB_1="LibB-220823-293T-1",
        libB_2="LibB-220823-293T-2",
        libB_3="LibB-220907-293T-1",
        libB_4="LibB-220907-293T-2",
        summary_dir="results/func_scores/",
        score_dir="results/func_scores/",
        html_dir="results/averaged_func_scores_ridgeplot/"
    output:
        html_output="results/averaged_func_scores_ridgeplot/averaged_func_scores_ridgeplot.html",
        nb="results/notebooks/visualize_func_scores_by_variant_type.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/averaged_ridgeplot_func_scores.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p libA_1 {params.libA_1} \
            -p libA_2 {params.libA_2} \
            -p libA_3 {params.libA_3} \
            -p libA_4 {params.libA_4} \
            -p libB_1 {params.libB_1} \
            -p libB_2 {params.libB_2} \
            -p libB_3 {params.libB_3} \
            -p libB_4 {params.libB_4} \
            -p summary_dir {params.summary_dir} \
            -p score_dir {params.score_dir} \
            -p html_dir {params.html_dir} \
            -p html_output {output.html_output} \
            &> {log}
        """


rule human_mastomys_correlation:
    """
    Correlation of functional scores from humanDAG1
    vs mastomysDAG1 expressing cells.
    """
    input:
        hek293T_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        humanDAG1_scores="results/func_effects/averages/human_293T_entry_func_effects.csv",
        mastomysDAG1_scores="results/func_effects/averages/mastomys_293T_entry_func_effects.csv",
        nb="notebooks/human_mastomys_correlation.ipynb",
    params:
        MTS=2, # min_times_seen
        n_selections=8,
        html_dir="results/DAG1_ortholog_correlations/"
    output:
        html_output="results/DAG1_ortholog_correlations/DAG1_ortholog_correlations.html",
        nb="results/notebooks/human_mastomys_correlation.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/human_mastomys_correlation.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p HEK293T_data_path {input.hek293T_scores} \
            -p humanDAG1_data_path {input.humanDAG1_scores} \
            -p mastomysDAG1_data_path {input.mastomysDAG1_scores} \
            -p MTS {params.MTS} \
            -p n_selections {params.n_selections} \
            -p html_dir {params.html_dir} \
            -p html_output {output.html_output} \
            &> {log}
        """

rule validation_titers:
    """
    Correlation of single mutant validations titers to 
    predicted values from DMS pipeline.
    """
    input:
        titers="data/single_mutant_functional_validations.csv",
        func_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        nb="notebooks/validation_titers.ipynb",
    params:
        out_dir="results/validation_plots/",
    output:
        saved_image_path="results/validation_plots/functional_validation_correlation.svg",
        nb="results/notebooks/validation_titers.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/validation_titers.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p validation_titers_path {input.titers} \
            -p functional_scores_path {input.func_scores} \
            -p out_dir {params.out_dir} \
            -p saved_image_path {output.saved_image_path} \
            &> {log}
        """

rule validation_neuts_89F:
    """
    Correlation of single mutant validations for 8.9F to 
    predicted values from DMS pipeline.
    """
    input:
        frac_infected="data/validation_frac_infected_89F.csv",
        neuts="data/validation_neuts_89F.csv",
        predicted_scores="results/antibody_escape/averages/89F_mut_effect.csv",
        nb="notebooks/validation_neuts_89F.ipynb",
    params:
        out_dir="results/validation_plots/",
    output:
        neuts_image_path="results/validation_plots/validation_neut_curves_89F.svg",
        corr_image_path="results/validation_plots/89F_validation_correlation.svg",
        nb="results/notebooks/validation_neuts_89F.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/validation_neuts_89F.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p fraction_infected_89F_path {input.frac_infected} \
            -p validation_neuts_89F_path {input.neuts} \
            -p model_predictions_path {input.predicted_scores} \
            -p out_dir {params.out_dir} \
            -p neuts_image_path {output.neuts_image_path} \
            -p corr_image_path {output.corr_image_path} \
            &> {log}
        """


rule visualize_RBD_regions:
    """
    Heatmaps of alpha-dystroglycan and LAMP1
    binding residues.
    """
    input:
        func_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        nb="notebooks/visualize_RBD_regions.ipynb",
    params:
        min_times_seen=2,
        n_selections=8,
        html_dir="results/func_scores_distributions/",
    output:
        html_output="results/func_scores_distributions/func_scores_distributions.html",
        nb="results/notebooks/visualize_RBD_regions.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/visualize_RBD_regions.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p func_scores {input.func_scores} \
            -p min_times_seen {params.min_times_seen} \
            -p n_selections {params.n_selections} \
            -p html_dir {params.html_dir} \
            -p html_output {output.html_output}
            &> {log}
        """


rule get_filtered_escape_CSVs:
    """
    Filter antibody escape data based on configuration 
    for displaying data to easily port to other analysis.
    """
    input:
        func_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        escape_377H="results/antibody_escape/averages/377H_mut_effect.csv",
        escape_89F="results/antibody_escape/averages/89F_mut_effect.csv",
        escape_2510C="results/antibody_escape/averages/2510C_mut_effect.csv",
        escape_121F="results/antibody_escape/averages/121F_mut_effect.csv",
        escape_256A="results/antibody_escape/averages/256A_mut_effect.csv",
        escape_372D="results/antibody_escape/averages/372D_mut_effect.csv",
        contacts_89F="data/antibody_contacts/antibody_contacts_89F.csv",
        contacts_377H="data/antibody_contacts/antibody_contacts_377H.csv",
        contacts_256A="data/antibody_contacts/antibody_contacts_256A.csv",
        contacts_2510C="data/antibody_contacts/antibody_contacts_2510C.csv",
        contacts_121F="data/antibody_contacts/antibody_contacts_121F.csv",
        contacts_372D="data/antibody_contacts/antibody_contacts_372D.csv",
        nb="notebooks/get_filtered_CSVs.ipynb",
    params:
        min_times_seen=2,
        min_func_score=-1.5,
        n_selections=8,
        frac_models=1,
        out_dir="results/filtered_antibody_escape_CSVs/",
        out_dir_images="results/antibody_escape_profiles/",
    output:
        filtered_escape_377H="results/filtered_antibody_escape_CSVs/377H_filtered_mut_effect.csv",
        filtered_escape_89F="results/filtered_antibody_escape_CSVs/89F_filtered_mut_effect.csv",
        filtered_escape_2510C="results/filtered_antibody_escape_CSVs/2510C_filtered_mut_effect.csv",
        filtered_escape_121F="results/filtered_antibody_escape_CSVs/121F_filtered_mut_effect.csv",
        filtered_escape_256A="results/filtered_antibody_escape_CSVs/256A_filtered_mut_effect.csv",
        filtered_escape_372D="results/filtered_antibody_escape_CSVs/372D_filtered_mut_effect.csv",
        func_effect_scale_bar="results/antibody_escape_profiles/func_effect_scale_bar.svg",
        escape_scale_bar="results/antibody_escape_profiles/escape_scale_bar.svg",
        saved_image_path="results/antibody_escape_profiles/antibody_escape_profiles.svg",
        validation_image_path="results/antibody_escape_profiles/validation_escape_profile.svg",
        nb="results/notebooks/get_filtered_CSVs.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/get_filtered_escape_CSVs.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p func_scores {input.func_scores} \
            -p escape_377H {input.escape_377H} \
            -p escape_89F {input.escape_89F} \
            -p escape_2510C {input.escape_2510C} \
            -p escape_121F {input.escape_121F} \
            -p escape_256A {input.escape_256A} \
            -p escape_372D {input.escape_372D} \
            -p contacts_89F {input.contacts_89F} \
            -p contacts_377H {input.contacts_377H} \
            -p contacts_256A {input.contacts_256A} \
            -p contacts_2510C {input.contacts_2510C} \
            -p contacts_121F {input.contacts_121F} \
            -p contacts_372D {input.contacts_372D} \
            -p min_times_seen {params.min_times_seen} \
            -p min_func_score {params.min_func_score} \
            -p n_selections {params.n_selections} \
            -p frac_models {params.frac_models} \
            -p out_dir {params.out_dir} \
            -p out_dir_images {params.out_dir_images} \
            -p filtered_escape_377H {output.filtered_escape_377H} \
            -p filtered_escape_89F {output.filtered_escape_89F} \
            -p filtered_escape_2510C {output.filtered_escape_2510C} \
            -p filtered_escape_121F {output.filtered_escape_121F} \
            -p filtered_escape_256A {output.filtered_escape_256A} \
            -p filtered_escape_372D {output.filtered_escape_372D} \
            -p func_effect_scale_bar {output.func_effect_scale_bar} \
            -p escape_scale_bar {output.escape_scale_bar} \
            -p saved_image_path {output.saved_image_path} \
            -p validation_image_path {output.validation_image_path} 
            &> {log}
        """

rule escape_sites_stratified_by_antibody_distance:
    """
    Stratify escape sites by distance to antibody.
    """
    input:
        contacts_89F="data/antibody_contacts/antibody_contacts_89F.csv",
        contacts_377H="data/antibody_contacts/antibody_contacts_377H.csv",
        contacts_256A="data/antibody_contacts/antibody_contacts_256A.csv",
        contacts_2510C="data/antibody_contacts/antibody_contacts_2510C.csv",
        contacts_121F="data/antibody_contacts/antibody_contacts_121F.csv",
        contacts_372D="data/antibody_contacts/antibody_contacts_372D.csv",
        filtered_escape_377H="results/filtered_antibody_escape_CSVs/377H_filtered_mut_effect.csv",
        filtered_escape_89F="results/filtered_antibody_escape_CSVs/89F_filtered_mut_effect.csv",
        filtered_escape_2510C="results/filtered_antibody_escape_CSVs/2510C_filtered_mut_effect.csv",
        filtered_escape_121F="results/filtered_antibody_escape_CSVs/121F_filtered_mut_effect.csv",
        filtered_escape_256A="results/filtered_antibody_escape_CSVs/256A_filtered_mut_effect.csv",
        filtered_escape_372D="results/filtered_antibody_escape_CSVs/372D_filtered_mut_effect.csv",
        func_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        nb="notebooks/escape_vs_antibody_distance.ipynb",
    params:
        out_dir="results/antibody_escape_profiles/",
    output:
        saved_image_path="results/antibody_escape_profiles/antibody_escape_by_distance.svg",
        nb="results/notebooks/escape_vs_antibody_distance.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/escape_sites_stratified_by_antibody_distance.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p contacts_89F {input.contacts_89F} \
            -p contacts_377H {input.contacts_377H} \
            -p contacts_256A {input.contacts_256A} \
            -p contacts_2510C {input.contacts_2510C} \
            -p contacts_121F {input.contacts_121F} \
            -p contacts_372D {input.contacts_372D} \
            -p filtered_escape_377H {input.filtered_escape_377H} \
            -p filtered_escape_89F {input.filtered_escape_89F} \
            -p filtered_escape_2510C {input.filtered_escape_2510C} \
            -p filtered_escape_121F {input.filtered_escape_121F} \
            -p filtered_escape_256A {input.filtered_escape_256A} \
            -p filtered_escape_372D {input.filtered_escape_372D} \
            -p func_scores {input.func_scores} \
            -p out_dir {params.out_dir} \
            -p saved_image_path {output.saved_image_path}
            &> {log}
        """


rule map_scores_onto_pdb_structure:
    """
    Map filtered functional and antibody scores onto pdb structure.
    """
    input:
        func_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        pdb_file="data/7puy.pdb",
        filtered_escape_377H="results/filtered_antibody_escape_CSVs/377H_filtered_mut_effect.csv",
        filtered_escape_89F="results/filtered_antibody_escape_CSVs/89F_filtered_mut_effect.csv",
        filtered_escape_2510C="results/filtered_antibody_escape_CSVs/2510C_filtered_mut_effect.csv",
        filtered_escape_121F="results/filtered_antibody_escape_CSVs/121F_filtered_mut_effect.csv",
        filtered_escape_256A="results/filtered_antibody_escape_CSVs/256A_filtered_mut_effect.csv",
        filtered_escape_372D="results/filtered_antibody_escape_CSVs/372D_filtered_mut_effect.csv",
        nb="notebooks/pdb_mapping.ipynb",
    params:
        min_times_seen=2,
        n_selections=8,
        out_dir="results/mapped_scores_onto_pdb/",
    output:
        pdb_func_min="results/mapped_scores_onto_pdb/func_scores_min.pdb",
        pdb_func_max="results/mapped_scores_onto_pdb/func_scores_max.pdb",
        pdb_377H="results/mapped_scores_onto_pdb/377H_escape.pdb",
        pdb_89F="results/mapped_scores_onto_pdb/89F_escape.pdb",
        pdb_2510C="results/mapped_scores_onto_pdb/2510C_escape.pdb",
        pdb_121F="results/mapped_scores_onto_pdb/121F_escape.pdb",
        pdb_256A="results/mapped_scores_onto_pdb/256A_escape.pdb",
        pdb_372D="results/mapped_scores_onto_pdb/372D_escape.pdb",
        nb="results/notebooks/pdb_mapping.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/map_scores_onto_pdb_structure.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p func_scores {input.func_scores} \
            -p pdb_file {input.pdb_file} \
            -p filtered_escape_377H {input.filtered_escape_377H} \
            -p filtered_escape_89F {input.filtered_escape_89F} \
            -p filtered_escape_2510C {input.filtered_escape_2510C} \
            -p filtered_escape_121F {input.filtered_escape_121F} \
            -p filtered_escape_256A {input.filtered_escape_256A} \
            -p filtered_escape_372D {input.filtered_escape_372D} \
            -p min_times_seen {params.min_times_seen} \
            -p n_selections {params.n_selections} \
            -p out_dir {params.out_dir} \
            -p pdb_func_min {output.pdb_func_min} \
            -p pdb_func_max {output.pdb_func_max} \
            -p pdb_377H {output.pdb_377H} \
            -p pdb_89F {output.pdb_89F} \
            -p pdb_2510C {output.pdb_2510C} \
            -p pdb_121F {output.pdb_121F} \
            -p pdb_256A {output.pdb_256A} \
            -p pdb_372D {output.pdb_372D} 
            &> {log}
        """


docs["Additional analyses and data files"] = {
    "Averaged functional score distributions by variant type" : {
        "Interactive plot showing averaged functional score distributions by variant type" : rules.averaged_ridgeplot_func_scores.output.html_output,
        "Notebook creating averaged functional score distributions" : rules.averaged_ridgeplot_func_scores.output.nb,
    },
    "Correlations of DAG1 ortholog functional selections" : {
        "Interactive plot showing correlations between DAG1 ortholog functional selections" : rules.human_mastomys_correlation.output.html_output,
        "Notebook correlating functional effects of humanDAG1 vs mastomysDAG1 expressing cells" : rules.human_mastomys_correlation.output.nb,
    },
    "Functional and antibody selection validations" : {
        "Notebook correlating measured vs predicted functional effects" : rules.validation_titers.output.nb,
        "Notebook correlating measured vs predicted neutralization": rules.validation_neuts_89F.output.nb,
    },
    "Functional scores for different GPC regions" : {
        "Interactive plot showing functional scores for different GPC regions" : rules.visualize_RBD_regions.output.html_output,
        "Notebook visualizing functional scores for different GPC regions" : rules.visualize_RBD_regions.output.nb,
    },
    "Filtered antibody escape data" : {
        "Notebook applying filters to antibody escape data" : rules.get_filtered_escape_CSVs.output.nb,
        "Filtered antibody escape CSVs" : {
            "CSV with filtered 377H antibody escape data" : rules.get_filtered_escape_CSVs.output.filtered_escape_377H,
            "CSV with filtered 89F antibody escape data" : rules.get_filtered_escape_CSVs.output.filtered_escape_89F,
            "CSV with filtered 2510C antibody escape data" : rules.get_filtered_escape_CSVs.output.filtered_escape_2510C,
            "CSV with filtered 121F antibody escape data" : rules.get_filtered_escape_CSVs.output.filtered_escape_121F,
            "CSV with filtered 256A antibody escape data" : rules.get_filtered_escape_CSVs.output.filtered_escape_256A,
            "CSV with filtered 372D antibody escape data" : rules.get_filtered_escape_CSVs.output.filtered_escape_372D,
        },
    },
    "Antibody escape stratified by distance to antibody" : {
        "Notebook plotting escape by distance to antibody" : rules.escape_sites_stratified_by_antibody_distance.output.nb,
    },
    "Mapped data onto pdb structure" : {
        "Notebook mapping score onto pdb structure" : rules.map_scores_onto_pdb_structure.output.nb,
    },
}
