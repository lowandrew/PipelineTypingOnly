#!/usr/bin/env python3
import subprocess
from accessoryFunctions.accessoryFunctions import *
from accessoryFunctions import metadataprinter
import glob

__author__ = 'adamkoziol', 'andrewlow'


def get_genome_info(fastafile):
    # TODO: Percent GC is giving a (slightly) different number than what qualimap was getting before. Look into this.
    # Looks to be because qualimap was giving GC content of mapped reads - this is probably a better way of doing things
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    num_bases = 0
    num_contigs = 0
    genome = SeqIO.parse(fastafile, 'fasta')
    total_seq = ""
    for contig in genome:
        num_contigs += 1
        num_bases += len(contig.seq)
        total_seq += str(contig.seq)
    gc_percent = str('%.2f' % GC(total_seq))
    return str(num_bases) + 'bp', str(num_contigs), gc_percent


def remove_unnecessary_columns(result_file):
    import csv
    import os
    necessary_columns = ['SampleName', 'N50', 'NumContigs', 'TotalLength', 'ReferenceGenome', 'RefGenomeAlleleMatches',
                         '16sPhylogeny', 'rMLSTsequenceType', 'MLSTsequencetype', 'MLSTmatches', 'coreGenome',
                         'Serotype', 'geneSeekrProfile', 'vtyperProfile', 'percentGC', 'TotalPredictedGenes',
                         'predictedgenesover3000bp', 'predictedgenesover1000bp', 'predictedgenesover500bp',
                         'predictedgenesunder500bp']
    with open(result_file.replace("combinedMetadata", "typingInfo"), 'w') as out_handle:
        tmp_row = list()
        for column in necessary_columns:
            tmp_row.append(column)
        out_handle.write(','.join(tmp_row) + "\n")
        with open(result_file, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                tmp_row = list()
                for column in necessary_columns:
                    tmp_row.append(row[column])
                out_handle.write(','.join(tmp_row) + "\n")
    os.remove(result_file)


class RunTyping(object):

    def quality(self):
        """
        Creates sequence objects and runs the quality assessment (FastQC), and quality trimming (bbduk) on the
        supplied sequences
        """
        from spadespipeline import mMLST
        from spadespipeline import quaster
        from spadespipeline import prodigal
        from spadespipeline import createobject
        import shutil
        # Create the objects
        self.runmetadata = createobject.ObjectCreation(self)
        # Determine the amount of physical memory in the system
        from psutil import virtual_memory
        mem = virtual_memory()
        # If the total amount of memory is greater than 100GB (this could probably be lowered), run CLARK
        if mem.total >= 100000000000:
            # Run CLARK typing on the .fastq and .fasta files
            from metagenomefilter import automateCLARK
            automateCLARK.PipelineInit(self, 'typing')
        else:
            printtime('Not enough RAM to run CLARK!', self.start)
        # Create a list of the files to process
        fasta_files = glob.glob(self.path + "*.fasta")
        # For each of the fasta files, create a folder based on the name of the file - 2014-SEQ-0276.fasta will
        # be placed in a folder named 2014-SEQ-0276
        for fasta in fasta_files:
            # Set the name of the folder
            fasta_dir = os.path.splitext(fasta)[0]
            # Create the folder
            make_path(fasta_dir)
            # Determine the name and extension of the fasta file
            fastaname = os.path.split(fasta)[-1]
            # Copy the file into the appropriate folder
            shutil.copy(fasta, os.path.join(fasta_dir, fastaname))
        # Perform quality assessments. After each method completes, print the metadata to file
        # Run gene predictions
        prodigal.Prodigal(self)
        metadataprinter.MetadataPrinter(self)
        # Run rMLST
        mMLST.PipelineInit(self, 'rmlst')
        metadataprinter.MetadataPrinter(self)
        # Run quast assembly metrics
        for sample in self.runmetadata.samples:
            sample.general.filteredfile = sample.general.bestassemblyfile
        quaster.Quast(self)
        # Calculate the depth of coverage as well as other quality metrics using Qualimap
        metadataprinter.MetadataPrinter(self)

    def typing(self):
        """
        Process samples with typing methods
        """
        from spadespipeline import mMLST
        from spadespipeline import GeneSeekr
        from spadespipeline import sixteenS
        from spadespipeline import univec
        from spadespipeline import prophages
        from spadespipeline import plasmidfinder
        from spadespipeline import serotype
        from spadespipeline import virulence
        from spadespipeline import vtyper
        from coreGenome import core
        from spadespipeline import sistr
        # Run modules and print metadata to file
        mMLST.PipelineInit(self, 'mlst')
        metadataprinter.MetadataPrinter(self)
        geneseekr = GeneSeekr.PipelineInit(self, 'geneseekr', True, 50, False)
        GeneSeekr.GeneSeekr(geneseekr)
        metadataprinter.MetadataPrinter(self)
        sixteens = GeneSeekr.PipelineInit(self, 'sixteenS', False, 95, False)
        sixteenS.SixteenS(sixteens)
        metadataprinter.MetadataPrinter(self)
        uni = univec.PipelineInit(self, 'univec', False, 80, False)
        univec.Univec(uni)
        metadataprinter.MetadataPrinter(self)
        pro = GeneSeekr.PipelineInit(self, 'prophages', False, 80, True)
        prophages.Prophages(pro)
        metadataprinter.MetadataPrinter(self)
        plasmid = GeneSeekr.PipelineInit(self, 'plasmidfinder', False, 80, True)
        plasmidfinder.PlasmidFinder(plasmid)
        metadataprinter.MetadataPrinter(self)
        sero = GeneSeekr.PipelineInit(self, 'serotype', True, 95, False)
        serotype.Serotype(sero)
        metadataprinter.MetadataPrinter(self)
        vir = GeneSeekr.PipelineInit(self, 'virulence', True, 70, True)
        virulence.Virulence(vir)
        # Remove samples that are not Escherichia so virulence finder doesn't attempt to work on them.
        to_be_readded = list()
        to_be_deleted = list()
        for sample in vir.runmetadata.samples:
            print(sample.general.referencegenus)
            if sample.general.referencegenus != 'Escherichia':
                to_be_readded.append(sample)
                to_be_deleted.append(sample)
        for sample in to_be_deleted:
            vir.runmetadata.samples.remove(sample)
        if len(vir.runmetadata.samples) > 0:
            GeneSeekr.GeneSeekr(vir)
        for sample in to_be_readded:
            self.runmetadata.samples.append(sample)
        metadataprinter.MetadataPrinter(self)
        # armiobject = GeneSeekr.PipelineInit(self, 'ARMI', False, 70)
        # armi.ARMI(armiobject)
        metadataprinter.MetadataPrinter(self)
        for sample in self.runmetadata.samples:
            for key in sample.geneseekr.blastresults:
                if "vt" in key or "VT" in key:
                    sample.general.stx = True
                    sample.general.filenoext = sample.general.filteredfile.split('.')[0]
        vtyper.Vtyper(self, 'vtyper')
        metadataprinter.MetadataPrinter(self)
        coregen = GeneSeekr.PipelineInit(self, 'coregenome', True, 70, False)
        core.CoreGenome(coregen)
        # Process core genome Escherichia samples
        core.AnnotatedCore(self)
        # metadataprinter.MetadataPrinter(self)
        sistr.Sistr(self, 'sistr')
        res = GeneSeekr.PipelineInit(self, 'resfinder', False, 80, True)
        GeneSeekr.GeneSeekr(res)
        # resfinder.ResFinder(res)
        metadataprinter.MetadataPrinter(self)

    def __init__(self, args, pipelinecommit, startingtime, scriptpath):
        """
        :param args: list of arguments passed to the script
        Initialises the variables required for this class
        """
        from spadespipeline import reporter
        from spadespipeline import compress
        import multiprocessing
        from spadespipeline import versions
        import time
        printtime('Welcome to the CFIA bacterial typing pipeline', self.starttime)
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.path, '')
        self.sequencepath = self.path
        self.start = time.time()
        self.reffilepath = os.path.join(args.referencefilepath, '')
        self.kmers = args.kmerrange
        self.updatedatabases = args.updatedatabases
        self.starttime = startingtime
        self.customsamplesheet = args.customsamplesheet
        if self.customsamplesheet:
            assert os.path.isfile(self.customsamplesheet), 'Cannot find custom sample sheet as specified {}'\
                .format(self.customsamplesheet)
        # Use the argument for the number of threads to use, or default to the number of cpus in the system
        self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
        # Assertions to ensure that the provided variables are valid
        make_path(self.path)
        assert os.path.isdir(self.path), u'Supplied path location is not a valid directory {0!r:s}'.format(self.path)
        self.reportpath = '{}reports'.format(self.path)
        assert os.path.isdir(self.reffilepath), u'Reference file path is not a valid directory {0!r:s}'\
            .format(self.reffilepath)
        self.commit = str(pipelinecommit)
        self.homepath = scriptpath
        self.runinfo = ''
        self.logfile = os.path.join(self.path, 'logfile.txt')
        # Initialise the metadata object
        self.runmetadata = MetadataObject()
        try:
            # Create the quality object
            self.quality()
            # Print the metadata to file
            metadataprinter.MetadataPrinter(self)
            # Perform typing of assemblies
            self.typing()
        except KeyboardInterrupt:
            raise KeyboardInterruptError
        # Create a report
        # Things reporter expects that I don't have here: NumContigs, TotalLength, MeanInsertSize, AverageCoverageDepth
        # All Run info. May need to modify the reporter.
        # This is an awful lot of dummy info. Some can be found manually (numcontigs and whatnot would be easy)
        # Things I want to actually get: Number of contigs, number of bases, gc percentage.
        for sample in self.runmetadata.samples:
            sample.mapping = GenObject()
            sample.mapping.Bases, sample.mapping.Contigs, sample.mapping.GcPercentage \
                = get_genome_info(sample.general.bestassemblyfile)
            # sample.mapping.Contigs = '1'
            # sample.mapping.Bases = "100000bp1"
            sample.mapping.MeanInsertSize = "200"
            sample.mapping.MeanCoveragedata = "30X"
            # sample.mapping.GcPercentage = "50"
            # sample.coregenome.targetspresent = '3'
            # sample.coregenome.totaltargets = '4'
            sample.run = GenObject()
            sample.run.Date = "2017-07-27"
            sample.run.InvestigatorName = "Darles Charwin"
            sample.run.TotalClustersinRun = '2'
            sample.run.NumberofClustersPF = '2'
            sample.run.PercentOfClusters = '1'
            sample.run.forwardlength = '1'
            sample.run.reverselength = '1'
            sample.run.SampleProject = '1'
            sample.general.spadesoutput = "dummyfolder"
        reporter.Reporter(self)
        compress.Compress(self)
        # Get all the versions of the software used
        versions.Versions(self)
        metadataprinter.MetadataPrinter(self)
        remove_unnecessary_columns(self.reportpath + "/combinedMetadata.csv")

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from time import time
    # from .accessoryFunctions import printtime
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git tag | tail -n 1'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Type bacterial strains from FASTA files.')
    parser.add_argument('-v', '--version',
                        action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',
                        help='Specify path')
    parser.add_argument('-t', '--threads',
                        help='Number of threads. Default is the number of cores in the system')
    parser.add_argument('-r', '--referencefilepath',
                        default='/spades_pipeline/SPAdesPipelineFiles',
                        help='Provide the location of the folder containing the pipeline accessory files (reference '
                             'genomes, MLST data, etc.')
    parser.add_argument('-k', '--kmerrange',
                        default='21,33,55,77,99,127',
                        help='The range of kmers used in SPAdes assembly. Default is 21,33,55,77,99,127')
    parser.add_argument('-c', '--customsamplesheet',
                        help='Path of folder containing a custom sample sheet and name of sample sheet file '
                             'e.g. /home/name/folder/BackupSampleSheet.csv. Note that this sheet must still have the '
                             'same format of Illumina SampleSheet.csv files')
    parser.add_argument('-u', '--updatedatabases',
                        action='store_true',
                        help='Optionally update (r)MLST databases')

    # Get the arguments into an object
    arguments = parser.parse_args()
    starttime = time()
    # Run the pipeline
    RunTyping(arguments, commit, starttime, homepath)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (starttime() - arguments.start) + '\033[0m')
