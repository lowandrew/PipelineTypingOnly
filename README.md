# PipelineTypingOnly

- This script runs the post-assembly parts of the OLC Pipeline. 

### Requirements
#### Programs you'll need installed:
- Quast (version 4.5)
- blast
- prodigal
- NCBI ePCR
- Python3

#### Python packages you'll need:
- see requirements.txt (use `pip3 install -r requirements.txt` to install all)

### Usage:
- To run: python3 Typing.py -r <reference_file_path> <assembly_folder>
- <reference_file_path> is where data files needed are stored. Currently at /mnt/nas/Adam/assemblypipeline/
- <assembly_folder> is the path to a folder containing assemblies to be typed, in FASTA format. To be safe, use absolute, not relative paths.


More to come...
