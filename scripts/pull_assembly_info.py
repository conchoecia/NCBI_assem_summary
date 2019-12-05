#!/usr/bin/env python

"""
This program takes a directory of assembly files from NCBI,
 gathers info from NCBI about the assembly, then runs The starting point for this script is this answer in biostars:
   https://www.biostars.org/p/345510/
"""
import argparse
from Bio import Entrez
from glob import glob
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import subprocess
import sys
import time

script_path = os.path.dirname(os.path.realpath(__file__))
fs_path = os.path.join(script_path, "../bin/fasta_stats")

field_list = ["GbUid", "SubmissionDate", "SpeciesName",
 "Organism", "FtpPath_GenBank", "SpeciesTaxid",
 "AssemblyAccession", "AssemblyType", "PartialGenomeRepresentation",
 "Coverage", "FtpPath_Assembly_rpt", "ChainId",
 "AssemblyClass", "AssemblyName", "FtpPath_Stats_rpt",
 "SortOrder", "AssemblyStatus", "BioSampleAccn",
 "Primary", "BioSampleId", "Taxid", "LastMajorReleaseAccession",
 "SubmitterOrganization"]
# If GB_BioProjects - keep everything
# if Synonym keep "Genbank"
SynList = ["Genbank", "AssemblyAccession"]

# just make sure that the path with the directories exists
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

#Arguments for the program
def argparser():
    """parses the arguments. Returns them as a dictionary."""
    parser = argparse.ArgumentParser(description='clustergenerator')
    parser.add_argument('-p', '--path',
                        type=dir_path,
                        required = True,
                        nargs='+',
                        help="""The directories that contains the assembly files.
                        Each directory name should be separated by a space.""")
    parser.add_argument('-e', '--email',
                        type=str,
                        required = True,
                        help="""The email associated with the NCBI account.""")
    parser.add_argument('-a', '--api_key',
                        type=str,
                        required = True,
                        help="""The api key associated with the email.""")
    args = parser.parse_args()
    return args

#Finds the ids associated with the assembly
def get_ids(term):
    ids = []
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    ids.append(record["IdList"])
    return ids

#Fetch raw output
def get_raw_assembly_summary(id):
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    #Return individual fields
    #XML output: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=79781&report=%22full%22
    #return(record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName']) #This will return the Assembly name
    return(record)

#JSON formatted output
def get_assembly_summary_json(id):
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    #Convert raw output to json
    return(json.dumps(record, sort_keys=True,indent=4, separators=(',', ': ')))

#dict formatted output
def get_assembly_summary_dict(id):
    this_dict = {}
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    temp = record["DocumentSummarySet"]["DocumentSummary"][0]
    for key in temp:
        if key in field_list:
            this_dict[key] = temp[key]
        elif key == "GB_BioProjects":
            if len(temp["GB_BioProjects"]) > 0:
                for subkey in temp["GB_BioProjects"][0]:
                    this_dict[subkey] = temp["GB_BioProjects"][0][subkey]
        elif key == "Synonym":
            for subkey in temp["Synonym"]:
                if subkey in SynList:
                    this_dict[subkey] = temp["Synonym"][subkey]
    #return dict output
    return(this_dict)

def path_to_filelist(path):
    """
    Gets a list of files to process based on the file endings.

    Acceptable file endings are:
      ["fna.gz", ".fna", ".fasta.gz", ".fasta", ".fa.gz", ".fa"]
    """
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    fts = ["fna.gz", ".fna", ".fasta.gz", ".fasta", ".fa.gz", ".fa"]
    final_files = []
    for tfile in onlyfiles:
        add_this=False
        for ending in fts:
            if tfile.endswith(ending):
                add_this = True
                break
        if add_this:
            final_files.append(os.path.join(path,tfile))
    return(final_files)

def run_fasta_stats(fpath):
    """
    This runs fasta_stats on an assembly and returns a dictionary of the values.
    """
    this_data = {}
    # This version is if you want N95_scaflen. Not a good idea for bad assemblies.
    #new_fields = ["num_scaffolds", "num_tigs", "tot_size_scaffolds",
    #              "tot_size_tigs", "scaffold_N50", "scaffold_L50",
    #              "contig_N50", "contig_L50", "perGap", "N95_scaflen"]
    new_fields = ["num_scaffolds", "num_tigs", "tot_size_scaffolds",
                  "tot_size_tigs", "scaffold_N50", "scaffold_L50",
                  "contig_N50", "contig_L50", "perGap"]
    tcmd = "{} {}".format(fs_path, fpath).split(" ")
    results = subprocess.run(tcmd, stdout=subprocess.PIPE).stdout.decode('utf-8').split()
    assert len(new_fields)==len(results)
    for i in range(len(results)):
        this_data[new_fields[i]] = results[i]
    return(this_data)

def main(args):
    #Increase query limit to 10/s & get warnings
    Entrez.email = args.email
    #Get one from https://www.ncbi.nlm.nih.gov/account/settings/ page
    Entrez.api_key= args.api_key
    #Get a file list of assemblies to process
    flist=[]
    if type(args.path) == str:
        # user has just passed one path
        flist = path_to_filelist(args.path)
    else:
        # user has passed multiple paths
        for tpath in args.path:
            templist = path_to_filelist(tpath)
            flist = flist + templist
    # for each file generate a line on the CSV
    all_samples = []
    for tfile in flist:
        print("Parsing: {}".format(tfile), file=sys.stderr)
        this_dict = {}
        this_dict["assem_file"] = tfile
        # get the NCBI info first
        term = "_".join(tfile.split("/")[-1].split("_")[:2])
        ids = get_ids(term)
        if len(ids[0]) == 0:
            for temp in ["AssemblyName", "Organism", "SpeciesName"]:
                this_dict[temp] = tfile.split("/")[-1]
        else:
            this_dict = get_assembly_summary_dict(ids[0][0])

        # now get the other assembly info from fasta_stats
        fstats_dict = run_fasta_stats(tfile)
        z = {**this_dict, **fstats_dict}
        all_samples.append(z)
        time.sleep(1)

    # now make a df with all the results
    df = pd.DataFrame(all_samples)
    df.to_csv("genome_info.tsv", sep='\t')

if __name__ == "__main__":
    args = argparser()
    sys.exit(main(args))
