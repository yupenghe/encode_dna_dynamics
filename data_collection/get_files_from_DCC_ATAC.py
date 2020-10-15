#!/usr/bin/env python2.7
import requests, json
import subprocess
import os

accessions = ["ENCSR012YAB",
              "ENCSR023QZX",
              "ENCSR031HDN",
              "ENCSR032HKE",
              "ENCSR068YGC",
              "ENCSR079GOY",
              "ENCSR088UYE",
              "ENCSR096JCC",
              "ENCSR102NGD",
              "ENCSR150EOO",
              "ENCSR150RMQ",
              "ENCSR154BXN",
              "ENCSR176BYZ",
              "ENCSR204ZTY",
              "ENCSR211OCS",
              "ENCSR217NOA",
              "ENCSR255XTC",
              "ENCSR261ICG",
              "ENCSR273UFV",
              "ENCSR282YTE",
              "ENCSR302LIV",
              "ENCSR310MLB",
              "ENCSR312LQX",
              "ENCSR335VJW",
              "ENCSR343TXK",
              "ENCSR358MOW",
              "ENCSR363SKQ",
              "ENCSR371KFW",
              "ENCSR377YDY",
              "ENCSR382RUC",
              "ENCSR384JBF",
              "ENCSR389CLN",
              "ENCSR451NAE",
              "ENCSR460BUL",
              "ENCSR465PYP",
              "ENCSR468GUI",
              "ENCSR486XAS",
              "ENCSR551WBK",
              "ENCSR552ABC",
              "ENCSR559FAJ",
              "ENCSR597BGP",
              "ENCSR603MWL",
              "ENCSR609OHJ",
              "ENCSR618HDK",
              "ENCSR623GSD",
              "ENCSR627OCR",
              "ENCSR652CNN",
              "ENCSR662KNY",
              "ENCSR668EIA",
              "ENCSR690VOH",
              "ENCSR700QBR",
              "ENCSR732OTZ",
              "ENCSR758IRM",
              "ENCSR785NEL",
              "ENCSR798FDL",
              "ENCSR810HQR",
              "ENCSR819QOJ",
              "ENCSR820ACB",
              "ENCSR836PUC",
              "ENCSR876SYO",
              "ENCSR896XIN",
              "ENCSR903GMO",
              "ENCSR961SMM",
              "ENCSR966ORC",
              "ENCSR976LWP",
              "ENCSR983JWA"
]

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

os.system("mkdir -p raw_ATAC/")

for accession in accessions:
    # accession = accession.rstrip()
    # This URL locates the ENCODE biosample with accession number ENCSR800JXR
    URL = "https://www.encodeproject.org/experiments/"+accession+"/?frame=embedded"

    # GET the object
    response = requests.get(URL, headers=HEADERS)
    
    # Extract the JSON response as a python dict
    related_dataset = response.json()

    # Get information
    assay_type = related_dataset['assay_term_name'].replace(" ","_")
    short_name = fullname2short[related_dataset['biosample_ontology']['term_name']]
    age = related_dataset['aliases'][0].split(":")[1].split("_")[0].upper()

    # print
    print("\t".join([age,short_name,assay_type,"",accession]))

    for ind_file in range(len(related_dataset['files'])):
        if related_dataset['files'][ind_file]['file_format'] != 'fastq':
            continue
        bio_rep = related_dataset['files'][ind_file]['biological_replicates']
        output_prefix = "_".join([age,short_name,str(bio_rep[0]),"ATAC"])
        url = 'https://www.encodeproject.org'+related_dataset['files'][ind_file]['href']
        output_filename = "raw_ATAC/" + output_prefix + "." + str(ind_file) + ".bam"
        subprocess.check_call(['curl',
                               '-RL',
                               url,
                               "-o",
                               output_filename
        ])
   
