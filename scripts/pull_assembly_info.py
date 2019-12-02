#!/usr/bin/python

"""
This program takes a directory of assembly files from NCBI,
 gathers info from NCBI about the assembly, then runs The starting point for this script is this answer in biostars:
   https://www.biostars.org/p/345510/
"""
import argparse, os
from Bio import Entrez

#Increase query limit to 10/s & get warnings
Entrez.email = 
#Get one from https://www.ncbi.nlm.nih.gov/account/settings/ page
Entrez.api_key=

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

term="GCA_003864495.1"

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
                        help="""The directory that contains the assembly files.""")
    args = parser.parse_args()
    args = vars(args)
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
    handle = Entrez.esummary(db="assembly",id=id,report="full")
    record = Entrez.read(handle)
    temp = record["DocumentSummarySet"]["DocumentSummary"][0]
    this_dict = {}
    for key in temp:
        print(key)
        if key in field_list:
            this_dict[key] = temp[key]
        elif key == "GB_BioProjects":
            for subkey in temp["GB_BioProjects"][0]:
                this_dict[subkey] = temp["GB_BioProjects"][0][subkey]
        elif key == "Synonym":
            for subkey in temp["Synonym"]:
                if subkey in SynList:
                    this_dict[subkey] = temp["Synonym"][subkey]
    #return dict output
    return(this_dict)
    
def main(args):
    #Test
    #for id in get_ids(term):
    #    this_rec = get_raw_assembly_summary(id) #For raw output
    #    this_summary = get_assembly_dict(this_rec)

if __name__ == "__main__":
    args = argparser()
    sys.exit(main(args))
