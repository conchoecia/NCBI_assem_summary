"""
This script downloads all of the refseq protozoa genomes with annotations
"""
import ftputil
from Bio import Entrez
from collections import defaultdict
from io import StringIO
import os
import urllib
import pandas as pd
import sys

def get_lineage(taxid, EMAIL, APIKEY):
    '''Get the full lineage for a taxid from Entrez'''
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    lineage = defaultdict(lambda: 'Unclassified')
    eai_to_lineage_dict = {2759:  "Eukaryota",
                           33090: "Plant",
                           33154: "Opisthokonta",
                           4751:  "Fungi",
                           28009: "Choanoflagellate",
                           33208: "Metazoa",
                           10197: "Ctenophora",
                           6040:  "Porifera",
                           6073:  "Cnidaria",
                           10226: "Placozoa",
                           33213: "Bilateria",
                           33511: "Deuterostomia",
                           33317: "Protostomia",
                           89593: "Craniata",
                           7711: "Chordata",
                           40674: "Mammal",
                           6231: "Nematoda",
                           88770: "Panarthropoda",
                           2697495: "Spiralia",
                           6447: "Mollusca",
                           1206794: "Ecdysozoa"}
    # set up all the lineage boolean values as False
    for key in eai_to_lineage_dict:
        lineage[eai_to_lineage_dict[key]] = False
    lineage["Organism"] = records[0]['ScientificName']

    for entry in records[0]['LineageEx']:
        if entry['Rank'] in ['superkingdom', 'kingdom', 'phylum',
                             'class', 'order', 'family',
                             'genus']:
            lineage[entry['Rank']] = entry['ScientificName'].split(' ')[-1]
            lineage["{}_id".format(entry['Rank'])] = int(entry['TaxId'].split(' ')[-1])
        # entry as integer
        eai = int(entry["TaxId"])
        if eai in eai_to_lineage_dict:
            lineage[eai_to_lineage_dict[eai]] = True

    # now figure out which group things belong to
    # Bilateria
    #   mammal
    #   non-mammal_craniata
    #   non-craniata_chordata
    #   non-chordata_deuterostome
    #   nematoda
    #   panarthropoda
    #   other_ecdysozoa
    #   mollusca
    #   other_spiralia
    #   other_ecdysozoa
    #   other_bilateria
    if lineage["Eukaryota"]:
        if lineage["Plant"]:
            lineage["group"] = "plants"
        elif lineage["Choanoflagellate"]:
            lineage["group"] = "choanoflagellate"
        elif lineage["Opisthokonta"]:
            if lineage["Fungi"]:
                lineage["group"] = "fungi"
            elif lineage["Choanoflagellate"]:
                lineage["group"] = "choanoflagellate"
            elif lineage["Metazoa"]:
                if lineage["Ctenophora"]:
                    lineage["group"] = "ctenophora"
                elif lineage["Porifera"]:
                    lineage["group"] = "porifera"
                elif lineage["Cnidaria"]:
                    lineage["group"] = "cnidaria"
                elif lineage["Placozoa"]:
                    lineage["group"] = "placozoa"
                elif lineage["Bilateria"]:
                    if lineage["Deuterostomia"]:
                        # first handle the deuterostomes
                        if lineage["Chordata"]:
                            if lineage["Craniata"]:
                                if lineage["Mammal"]:
                                    lineage["group"] = "mammal"
                                    pass
                                else:
                                    lineage["group"] = "non-mammal_craniata"
                            else:
                                lineage["group"] = "non-craniata_chordata"
                        else:
                            lineage["group"] = "non-chordata_deuterostome"
                    elif lineage["Protostomia"]:
                        # handle protostomes
                        if lineage["Ecdysozoa"]:
                            if lineage["Nematoda"]:
                                lineage["group"] = "nematoda"
                            elif lineage["Panarthropoda"]:
                                lineage["group"] = "panarthropoda"
                            else:
                                lineage["group"] = "other_ecdysozoa"
                        elif lineage["Spiralia"]:
                            if lineage["Mollusca"]:
                                lineage["group"] = "mollusca"
                            else:
                                lineage["group"] = "other_spiralia"
                        else:
                            lineage["group"] = "other_ecdysozoa"
                    else:
                        lineage["group"] = "other_bilateria"
            else:
                lineage["group"] = "HTC"
        else:
            lineage["group"] = "protist"
    else:
        raise IOError("Shouldn't have arrived here")
    return dict(lineage)
