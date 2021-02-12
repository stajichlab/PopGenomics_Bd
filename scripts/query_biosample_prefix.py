#!/usr/bin/env python3

import csv, re, sys, os
import xml.etree.ElementTree as ET
from Bio import Entrez
Entrez.email = 'jason.stajich@ucr.edu'
insampleslst = ["samples.csv", "samples_single.csv"]
outsamples="samples_metadata.csv"

if len(sys.argv) > 1:
    insamples = sys.argv[1]

if len(sys.argv) > 2:
    outsamples = sys.argv[2]

seen = {}
if os.path.exists(outsamples):
    with open(outsamples,"rU") as preprocess:
        incsv = csv.reader(preprocess,delimiter=",")
        h = next(incsv)
        for row in incsv:
            name = "%s"%(row[0])
            print('storing %s for previous see'%(name))
            seen[name] = row

with open(outsamples,"w") as outfh:
    for insamples in insampleslst:
        with open(insamples,"rU") as infh:
            outcsv    = csv.writer(outfh,delimiter=",")
            outcsv.writerow(['SPECIES','STRAIN','BIOSAMPLE','BIOPROJECT','SRA','LOCUSTAG',
                             'LOCATION','YEAR','HOST','LAT_LONG'])

            samplescsv = csv.reader(infh,delimiter=",")
            header = next(samplescsv)
            for row in samplescsv:
                strain = row[0]
                lookup = "%s"%(strain)
                outrow = [strain]
                print("strain is '%s' lname=%s"%(strain,lookup))
                if lookup in seen and len(seen[lookup][2]) > 0:
                    outrow = seen[lookup]
                    outcsv.writerow(outrow)
                    outfh.flush()
                    continue
                else:
                    print("doing a lookup for %s"%(lookup))

                handle = Entrez.esearch(db="biosample",retmax=10,term=lookup)
                record = Entrez.read(handle)
                handle.close()
                SRA = ""
                BIOSAMPLE = ""
                BIOPROJECT = ""
                LOCUSTAG = ""
                BIOPROJECTID=""
                LOCATION = ""
                LAT_LONG = ""
                for biosampleid in record["IdList"]:
                    handle = Entrez.efetch(db="biosample", id=biosampleid)
                    tree = ET.parse(handle)
                    root = tree.getroot()
                    for sample in root:
                        #print( ET.tostring(sample))
                        BIOSAMPLE = sample.attrib['accession']
                        for ids in root.iter('Ids'):
                            for id in ids.iter('Id'):
                                if 'db' in id.attrib and id.attrib['db'] == "SRA":
                                    SRA = id.text
                                elif 'db' not in id.attrib:
                                    print("missing db or SRA")
                                    #print( ET.tostring(root))
                                    for k in id.attrib:
                                        print("  have [%s] = '%s'"%(k,id.attrib[k]))

                        for links in root.iter('Links'):
                            for link in links:
                                linkdat = link.attrib
                                if 'type' in linkdat and linkdat['type'] == 'entrez' and 'label' in linkdat:
                                    BIOPROJECT = linkdat['label']
                                    BIOPROJECTID = link.text
                        for att in root.iter('Attributes'):
                            for a in att:
                                #print(a.attrib)
                                attname = a.attrib["attribute_name"]
                                #print(attname)
                                if 'lat_long' == attname:
                                    print(a.text)
                                    LAT_LONG = a.text
                                if 'geo_loc_name' == attname:
                                    LOCATION = a.text.encode('utf-8').strip()
                                    print(LOCATION)
                if BIOPROJECTID:
                    bioproject_handle = Entrez.efetch(db="bioproject",id = BIOPROJECTID)
                    projtree = ET.parse(bioproject_handle)
                    projroot = projtree.getroot()

                    lt = projroot.iter('LocusTagPrefix')
                    for locus in lt:
                        LOCUSTAG = locus.text
                outrow.extend([BIOSAMPLE,BIOPROJECT,SRA,LOCUSTAG,LOCATION,LAT_LONG])
                outcsv.writerow(outrow)
                outfh.flush()
