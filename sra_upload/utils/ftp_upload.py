import datetime
import ftplib
import os
import tarfile
import yaml
import argparse
import pandas as pd

def make_tar_file(fastq_file, filename):
    """
    Make the compressed *.tar file with runs to upload to SRA.

    Parameters
    ----------
    fastq_file: str
        Path to table of the fastq file names and paths for upload to the SRA.
    """
    
    # Reformat the fastqs so that there is one row per filename
    fastqs_long = pd.read_csv(fastq_file)
    fastqs = fastqs_long.groupby('filename', as_index=False).agg(" ".join)

    # Add all of the fastq files to upload to the tar directory
    try:
        
        # Make a temp directory to store concatenated fastqs
        if os.path.exists("tmp/"):
            os.system("rm -rf tmp/")
        os.makedirs("tmp/")

        # Open up the tar file
        with tarfile.open(filename, mode='w') as f:
            # Iterate over the unique filenames which can have many files per sample
            for i, tup in enumerate(fastqs.itertuples()):

                # For concise output, only alert after every 10 processed files. 
                if (i+1) % 10 == 0 or i == 0 or i+1 == len(fastqs):
                    print(f"Concatenating and Adding file {i + 1} of {len(fastqs)} to {filename}")

                # Concatenate multiple files into a single file corresponding to a library 
                os.system(f"cat {tup.filename_fullpath} > tmp/{tup.filename}")
                # Add the concatenated file to the tar
                f.add(f"tmp/{tup.filename}", arcname=tup.filename)

            print(f"Added all files to {filename}\n")
            
            # When done, remove the local copies of the fastq files.
            print("Removing tmp directory.\n")
            os.system("rm -rf tmp/")

    except:
        if os.path.isfile(filename):
           raise 

    print(f"The size of {filename} is {os.path.getsize(filename) / 1e9:.1f} GB\n")

    # Check that all of the files are in the tarred directory
    with tarfile.open(filename) as f:
        files_in_tar = set(f.getnames())
    if files_in_tar == set(fastqs['filename']):
        print(f"{filename} contains all {len(files_in_tar)} expected files.\n")
    else:
        raise Exception(f"{filename} does not have all the expected files.")
        
    print("Finished preparing the tar file for upload.")


def upload_via_ftp(tar_path,
                   ftp_username,
                   ftp_account_folder,
                   ftp_subfolder,
                   ftp_address='ftp-private.ncbi.nlm.nih.gov',
                   ftp_password='ftp_password.txt'
):
    """
    Upload the *.tar file of fastqs to the SRA with FTP.

    Parameters
    ----------
    tar_path: str
        Path to the tarred directory of fastqs.
    ftp_username: str
        Username provided by the SRA.
    ftp_account_folder: str
        Account folder provided by the SRA. 
    ftp_subfolder: str
        Unique sub-folder name chosen for upload.
    """
    
    # Get the password 
    with open(ftp_password) as f:
        password = f.read().strip()
        
    # Start the upload - this can take awhile
    print(f"Starting upload at {datetime.datetime.now()}")

    with ftplib.FTP(ftp_address) as ftp:
        ftp.login(user=ftp_username,
                  passwd=password,
                  )
        ftp.cwd(ftp_account_folder)
        ftp.mkd(ftp_subfolder)
        ftp.cwd(ftp_subfolder)
        with open(tar_path, 'rb') as f:
            ftp.storbinary(f"STOR {tar_path}", f)

    print(f"Finished upload at {datetime.datetime.now()}")


if __name__ == '__main__':
    """
    Have the option to run a command line version of the FTP upload because this can take 
    a really, really long time. 
    """

    # Command line interface
    parser = argparse.ArgumentParser(
        description="Submit samples to the SRA using the FTP protocol."
    )

    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the config file with FTP credentials.",
    )

    parser.add_argument(
        "--sampletype",
        type=str,
        required=True,
        help="Either barcodes or pacbio depending on what you're submitting.",
    )

    args = parser.parse_args()

    with open(args.config) as f:
         config = yaml.safe_load(f)
    
    print("Starting to upload samples to the SRA via FTP.\n")

    print(f"""
    ftp_username: {config['ftp_username']}
    ftp_account_folder: {config['ftp_account_folder']}
    """)

    if os.path.isfile('ftp_password.txt'):
        print("ftp_password.txt Exists!")
    else:
        raise Exception("Make sure that ftp_password.txt exists in this directory.")

    assert args.sampletype in {'pacbio', 'barcodes'}

    if args.sampletype == 'barcodes': 

        print("\nUploading runs for Illumina Barcodes.\n")
        
        print(f"ftp_subfolder: {config['barcode_runs']['ftp_subfolder']}")

        upload_via_ftp(tar_path="barcode_SRA_submission.tar",
        ftp_username=config['ftp_username'],
        ftp_account_folder=config['ftp_account_folder'],
        ftp_subfolder=config['barcode_runs']['ftp_subfolder'],
        ftp_address='ftp-private.ncbi.nlm.nih.gov',
        ftp_password='ftp_password.txt'
        )

    else: 

        print("\nUploading runs for PacBio CCSs.\n")
        
        print(f"ftp_subfolder: {config['pacbio_runs']['ftp_subfolder']}")

        upload_via_ftp(tar_path="pacbio_SRA_submission.tar",
        ftp_username=config['ftp_username'],
        ftp_account_folder=config['ftp_account_folder'],
        ftp_subfolder=config['pacbio_runs']['ftp_subfolder'],
        ftp_address='ftp-private.ncbi.nlm.nih.gov',
        ftp_password='ftp_password.txt'
        )