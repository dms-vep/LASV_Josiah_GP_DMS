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

rule human_mastomys_correlation:
    """
    Correlation of functional scores from humanDAG1
    vs mastomysDAG1 expressing cells.
    """
    input:
        humanDAG1_scores="results/func_effects/averages/human_293T_entry_func_effects.csv",
        mastomysDAG1_scores="results/func_effects/averages/mastomys_293T_entry_func_effects.csv",
        nb="notebooks/human_mastomys_correlation.ipynb",
    output:
        nb="results/notebooks/human_mastomys_correlation.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/human_mastomys_correlation.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            &> {log}
        """

rule validation_titers:
    """
    Correlation of single mutant validations titers to 
    predicted values from DMS pipeline.
    """
    input:
        titers="data/single_mutant_functional_validations.csv",
        predicted_scores="results/func_effects/averages/293T_entry_func_effects.csv",
        nb="notebooks/validation_titers.ipynb",
    output:
        nb="results/notebooks/validation_titers.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/validation_titers.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            &> {log}
        """

rule validation_neuts_89F:
    """
    Correlation of single mutant validations for 8.9F to 
    predicted values from DMS pipeline.
    """
    input:
        neuts="data/validation_neuts_89F.csv",
        predicted_scores="results/antibody_escape/averages/89F_mut_effect.csv",
        nb="notebooks/validation_neuts_89F.ipynb",
    output:
        nb="results/notebooks/validation_neuts_89F.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/validation_neuts_89F.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            &> {log}
        """


docs["Additional analyses"] = {
    "Additional analysis notebooks" : {
        "Notebook correlating functional effects on humanDAG1 vs mastomysDAG1 expressing cells" : rules.human_mastomys_correlation.output.nb,
        "Notebook correlating measured vs predicted functional effects" : rules.validation_titers.output.nb,
        "Notebook correlating measured vs predicted neutralization": rules.validation_neuts_89F.output.nb,
    },
}
