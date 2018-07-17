#!/usr/bin/env python2.7
import requests, json
import subprocess

fullname2short = {"embryonic facial prominence":"CF",
		  "forebrain":"FB",
		  "heart":"HT",
		  "hindbrain":"HB",
		  "intestine":"IT",
		  "kidney":"KD",
		  "limb":"LM",
		  "liver":"LV",
		  "lung":"LG",
		  "midbrain":"MB",
		  "neural tube":"NT",
		  "stomach":"ST"}

download = False
download = True

# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

for accession in file("atac.tsv").readlines():
    accession = accession.rstrip()
    # This URL locates the ENCODE biosample with accession number ENCSR800JXR
    URL = "https://www.encodeproject.org/experiments/"+accession+"/?frame=embedded"

    # GET the object
    response = requests.get(URL, headers=HEADERS)
    
    # Extract the JSON response as a python dict
    related_dataset = response.json()

    # Get information
    assay_type = related_dataset['assay_term_name'].replace(" ","_")
    short_name = fullname2short[related_dataset['biosample_term_name']]
    age = related_dataset['aliases'][0].split(":")[1].split("_")[0].upper()

    # print
    print("\t".join([age,short_name,assay_type,"",accession]))
