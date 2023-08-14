"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


rule validation_neuts_89F:
    """
    Correlation of single mutant validations for 8.9F to 
    predicted values from DMS pipeline.
    """
    input:
        neuts="data/validation_neuts_89F.csv",
        predicted_scores="results/antibody_escape/averages/89F_mut_icXX.csv",
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


docs["Single mutant validations"] = {
    "Notebook correlating measured vs predicted neutralization": rules.validation_neuts_89F.output.nb,
}
