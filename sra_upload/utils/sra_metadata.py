import os
import itertools 
import yaml
import pandas as pd


def format_metadata_tsv(config, outpath, pacbio=False):
    """ Format the SRA metadata tsv for a new BioSample.
    
    Parameters
    ----------
    config: dict
        Dict from a YAML format config file
    outpath: str
        Path to where outputs should be located
    pacbio: bool
        True if this BioSample is for PacBio reads

    """

    if not pacbio:
        
        # Extracting only the information for the illumina runs
        illumina_config = config['barcode_runs']
        illumina_runs = pd.read_csv(illumina_config['file_path'])
        
        # Check that the columns to make the library ID exist in the illumina runs table
        if not all(column in illumina_runs.columns for column in illumina_config['sample_id_columns']):
            raise Exception("Make sure that the columns specified in `sample_id_columns` exist in your barcoded runs table.")
        
        # Change NaNs to 'none' in any columns used for the library id name
        illumina_runs[illumina_config['sample_id_columns']] = illumina_runs[illumina_config['sample_id_columns']].fillna('none')
        
        # Make a new column for the sample name
        illumina_runs['sample'] = (illumina_runs[illumina_config['sample_id_columns']].astype(str) + '_').sum(axis=1).str.strip("_")
        
        # Add the necessary columns for submission 
        submission_metadata = (
            illumina_runs
            .assign(
                biosample_accession=lambda x: illumina_config['accession'],
                library_ID=lambda x: x['library'] + "_" + x['sample'],
                title=lambda x: illumina_config['title_prefix'] + " " + x['sample'],
                library_strategy=illumina_config['strategy'],
                library_source=illumina_config['source'],
                library_selection=illumina_config['selection'],
                library_layout=illumina_config['layout'],
                platform=illumina_config['platform'],
                instrument_model=illumina_config['model'],
                design_description=illumina_config['description'],
                filetype='fastq',
                filename_fullpath=lambda x: x['fastq_R1'].str.split(';')
            )
            .explode('filename_fullpath')
            .assign(
                filename=lambda x: x['library_ID']  + ".fastq.gz",
                filename_fullpath=lambda x: x['filename_fullpath'].str.strip()
            )
            .drop(columns=illumina_runs.columns)
            .reset_index(drop=True)
        )
        
        # All of the files should exist at their specified path
        files_for_upload = submission_metadata[['filename_fullpath', 'filename']]

        if not files_for_upload['filename_fullpath'].map(os.path.isfile).all():
            raise Exception("Check that the absolute filepaths to your sequencing runs are correct and exist.")
        
        # Widen the metadata so there are unqiue columns for each file associated with a single sample
        submission_wide_metadata = (
            submission_metadata
            .drop(columns=['filename_fullpath'])
            .drop_duplicates()
            )

        assert len(illumina_runs) == len(submission_wide_metadata)
        
    else:
        
        # Extracting only the information for the pacbio runs
        pacbio_config = config['pacbio_runs']
        pacbio_runs = pd.read_csv(pacbio_config['file_path'])
        
        # Add the necessary columns for submission 
        submission_metadata = (
            pacbio_runs
            .assign(
                biosample_accession=lambda x: pacbio_config['accession'],
                library_ID=lambda x: x['library'] + "_" + x['run'].astype(str) + "_PacBio_CCSs",
                title=lambda x: pacbio_config['title_prefix'] + " " + x['library_ID'],
                library_strategy=pacbio_config['strategy'],
                library_source=pacbio_config['source'],
                library_selection=pacbio_config['selection'],
                library_layout=pacbio_config['layout'],
                platform=pacbio_config['platform'],
                instrument_model=pacbio_config['model'],
                design_description=pacbio_config['description'],
                filetype='fastq',
                filename_fullpath=lambda x: x['fastq'],
                filename=lambda x: x['library_ID'] + '.fastq.gz'
            )
            .drop(columns=pacbio_runs.columns)
            .reset_index(drop=True)
        )
        
        # All of the files should exist at their specified path
        files_for_upload = submission_metadata[['filename_fullpath', 'filename']]

        if not files_for_upload['filename_fullpath'].map(os.path.isfile).all():
            raise Exception("Check that the absolute filepaths to your sequencing runs are correct and exist.")
        
        # All we need to do to finish formatting the PacBio is drop the file paths
        submission_wide_metadata = (
            submission_metadata
            .drop(columns=['filename_fullpath'])
            .drop_duplicates()
        )
        
    # Write out the metadata and corresponding files to be uploaded
    library_type = "barcodes" if not pacbio else "pacbio"
    submission_wide_metadata.to_csv(os.path.join(outpath, f"{library_type}_SRA_metadata.tsv"), sep='\t', index=False)
    files_for_upload.to_csv(os.path.join(outpath, f"{library_type}_fasta_files.csv"), index=False)


