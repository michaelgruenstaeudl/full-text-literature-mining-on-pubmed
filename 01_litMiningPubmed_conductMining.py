#!/usr/bin/env python3
__version__ = 'mi.gruenstaeudl@gmail.com|2022-08-04T22:09:36 CEST'

#-----------------------------------------------------------------#
## IMPORTS
import argparse
import Bio
from Bio import Entrez  # line is necessary, see: https://www.biostars.org/p/13099/
import bs4
import coloredlogs
import collections
import copy
import datetime
import ipdb
import json
import logging
import lxml
from lxml import etree  # line is necessary, see: https://stackoverflow.com/questions/41066480/lxml-error-on-windows-attributeerror-module-lxml-has-no-attribute-etree
import pprint
import re
import time
import urllib


#-----------------------------------------------------------------#
# INFO
#__info__ = 'Full text literature mining on PubMed'

#-----------------------------------------------------------------#
# DEBUGGING HELP
pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(my_dict)
#print(json.dumps(masterDict, indent=4))

#ipdb.set_trace()

#-----------------------------------------------------------------#
# CLASSES AND FUNCTIONS
class PubMedInteract:

    def __init__(self, email):
        self.email = email

    def search_via_query(self, query):
        '''Search PubMed via a query'''
        Bio.Entrez.email = self.email
        handle = Bio.Entrez.esearch(db='pubmed', 
                                sort='relevance', 
                                retmax='100000',  # set to 20 by default
                                #retmode='xml',  # set to xml by default
                                term=query)
        result = Bio.Entrez.read(handle)
        return result

    def retrieve_metadata(self, uid):
        '''Fetch metadata on PubMed for a given uid'''
        Bio.Entrez.email = self.email
        handle = Bio.Entrez.efetch(db='pubmed',  # Why does "db='pmc'" not work?
                               retmode='xml',
                               id=uid)
        result = Bio.Entrez.read(handle)
        return result

    def PMCID_to_PMID(self, pmid):
        '''Look up PubMedCentral ID from a PubMed ID'''
        Bio.Entrez.email = self.email
        handle = Bio.Entrez.elink(dbfrom='pubmed',
                              db='pmc',
                              linkname='pubmed_pmc',
                              id=pmid,
                              retmode='text')
        result = Bio.Entrez.read(handle)
        try:
            pmcid = 'PMC' + result[0]['LinkSetDb'][0]['Link'][0]['Id']
        except:
            pmcid = None
        return pmcid

class XMLparsing:

    def __init__(self):
        pass

    def article_title(self, inputxml):
        '''Parse the article title'''
        try:
            article_title = inputxml['MedlineCitation']['Article']['ArticleTitle']
        except:
            raise Exception("ERROR: Cannot parse the obligatory article title.")
        return article_title

    def article_publdate(self, inputxml):
        '''Parse the article year'''
        try:
            article_publdate = inputxml['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
        except:
            try:
                article_publdate = inputxml['MedlineCitation']['Article']['ArticleDate'][0]
            except:
                article_publdate = None
        return article_publdate

    def article_keywords(self, inputxml):
        '''Parse the article keywords'''
        try:
            #article_keywords = ' - '.join([p[:] for p in inputxml['MedlineCitation']['KeywordList'][0]])
            article_keywords = [p[:] for p in inputxml['MedlineCitation']['KeywordList'][0]]
        except:
            try:
                #article_keywords = ' - '.join([p['DescriptorName'][:] for p in inputxml['MedlineCitation']])
                article_keywords = [p['DescriptorName'][:] for p in inputxml['MedlineCitation']]
            except:
                article_keywords = None
        return article_keywords

    def article_mesh(self, inputxml):
        '''Parse the article MeSH (medical subject headings)'''
        try:
            #article_mesh = ' - '.join([p['DescriptorName'][:] for p in inputxml['MedlineCitation']['MeshHeadingList']])
            article_mesh = [p['DescriptorName'][:] for p in inputxml['MedlineCitation']['MeshHeadingList']]
        except:
            article_mesh = None
        return article_mesh

    def article_pmid(self, inputxml):
        '''Parse the article PMID (PubMed ID)'''
        try:
            article_pmid = ''.join([p[:] for p in inputxml['MedlineCitation']['PMID']])
        except:
            raise Exception("ERROR: Cannot parse the obligatory article PubMed ID (PMID).")
        return article_pmid

class PMCIDops:

    def __init__(self):
        pass
    
    def article_pmcid(self, article_pmid, email):
        '''Parse the article PMCID (PubMedCentral ID) from the PMID (PubMed ID)'''
        try:
            article_pmcid = PubMedInteract(email).PMCID_to_PMID(article_pmid)
        except:
            article_pmcid = None
        return article_pmcid
    
    def get_pmcid_bioc_url(self, article_pcmid):
        '''Form URL to access the fulltext via the PubMedCentral BioC API'''
        BioC_PMCID_pre = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/'
        BioC_PMCID_post = '/unicode/'  # Alternative: '/ascii/'
        return BioC_PMCID_pre + article_pcmid + BioC_PMCID_post

    def get_pmcid_pubreader_url(self, article_pcmid):
        '''Form URL to access the fulltext via the PubMedCentral PubReader'''
        PubReader_PMCID_pre = 'https://www.ncbi.nlm.nih.gov/pmc/articles/'
        PubReader_PMCID_post = '/?report=reader'
        return PubReader_PMCID_pre + article_pcmid + PubReader_PMCID_post

class LXMLops:

    def __init__(self, etree):
        self.etree = etree

    def remove_expendable(self):
        '''Remove unnecessary XML sections from full text'''
        xml_document = self.etree.find('.//document')
        # STEP 1. Removing figures, tables and any backmatter parts
        for passage in self.etree.findall('.//passage'):
            if passage.find('infon[@key="section_type"]').text in ['FIG', 'TABLE', 'COMP_INT', 'AUTH_CONT', 'SUPPL', 'ACK_FUND']:
                xml_document.remove(passage)  
        # STEP 2. Removing references
        for passage in self.etree.findall('.//passage'):
            if passage.find('infon[@key="type"]').text == 'ref':
                xml_document.remove(passage)

    def extract_all_text(self):
        '''Extract all text of the full text by paragraph'''
        paragraphs = []
        for passage in self.etree.findall('.//passage'):
            header = passage.find('infon[@key="section_type"]').text
            if header not in paragraphs:
                paragraphs.append(header)
            main_text = passage.find('text').text
            if main_text:
                paragraphs.append(main_text)
        #for paragr in paragraphs:
        #    for paragraph in re.findall(regex_for_sentence_delin, paragr):
        #        paragraphs.append(paragraph)
        return paragraphs

class HTMLops:

    def __init__(self):
        pass
    
    def remove_expendable(self, soup):
        '''Remove unnecessary sections from full text'''
        html_article = soup.find('article', {'data-type' : 'main'})
        if html_article.find('div', {'id' : 'ack-1'}):
            html_article.find('div', {'id' : 'ack-1'}).decompose()
        if html_article.find('div', {'id' : 'fn-group-1'}):
            html_article.find('div', {'id' : 'fn-group-1'}).decompose()
        if html_article.find('div', {'id' : '__ffn_sec'}):
            html_article.find('div', {'id' : '__ffn_sec'}).decompose()
        if html_article.find('div', {'id' : 'ref-list-1'}):
            html_article.find('div', {'id' : 'ref-list-1'}).decompose() 
        html_divs = html_article.find_all('div', {'class': 'tsec sec'})
        return html_divs

    def extract_all_text(self, in_list):
        '''Extract all text of the full text by paragraph'''
        all_paragraphs = []
        for div in in_list:
            all_paragraphs += [str(paragr.text) for paragr in div.find_all('p')]
        all_paragraphs = list(filter(None, all_paragraphs))  # removing empty strings in list
        return all_paragraphs

class ParagrOps:

    def __init__(self):
        pass
    
    def extract_gene_paragraphs(self, in_list, keywVariants_list):
        '''Extract only those paragraphs that contain the keyword `gene` in its different forms (and are, thus, potential hits)'''
        paragraphs = []
        for paragr in in_list:
            if any(word in paragr for word in keywVariants_list):
                paragraphs.append(paragr.strip())
        return paragraphs

    def remove_duplicate_paragraphs(self, in_list):
        '''Remove all duplicate paragraphs'''
        unique_paragraphs = []
        known_paragraphs = []
        for paragr in in_list:
            if paragr not in known_paragraphs:
                unique_paragraphs.append(paragr)
                known_paragraphs.append(paragr)
            else:
                pass
                #print("Duplicate: \n\n %s" % paragr)
        return unique_paragraphs

class MiscOps:

    def __init__(self):
        pass
    
    def convert_to_date(self, in_dict):
        '''Convert the date dictionary into a string that can be ordered'''
        outdate = []
        try:
            outdate.append(in_dict['Year'])
        except:
            outdate.append('1900')
        try:
            month = time.strptime(in_dict['Month'], '%b').tm_mon
            outdate.append(str(month).zfill(2))
        except:
            try:
                outdate.append(in_dict['Month'])
            except:
                outdate.append('01')
        return '-'.join(outdate)
    
    def write_to_json(self, json_filename, json_data):
        with open(json_filename, "w") as json_file:
            json.dump(json_data, json_file)
    

def main(args):

    ### STEP 1. Getting variables
    email = args.mail
    query = args.query
    verbose = args.verbose
    keywVariants_list = args.keywVariants_list
    keywVariants_list = [keyw.strip() for keyw in keywVariants_list]

    toxinDB_fname = args.toxinDB_fname
    with open(toxinDB_fname) as file:
        toxinDB_list = file.read().splitlines()
        toxinDB_list = [line.strip() for line in toxinDB_list]

    geneProdDB_fname = args.geneProdDB_fname
    with open(geneProdDB_fname) as file:
        geneProdDB_list = file.read().splitlines()
        geneProdDB_list = [line.strip() for line in geneProdDB_list]

    geneProdExcl_fname = args.geneProdExcl_fname
    with open(geneProdExcl_fname) as file:
        geneProdExcl_list = file.read().splitlines()
        geneProdExcl_list = [line.strip() for line in geneProdExcl_list]
        geneProdExcl_list = [line for line in geneProdExcl_list if '#' not in line]
        geneProdExcl_list = [line for line in geneProdExcl_list if not line.startswith('#')]
        geneProdExcl_list = list(filter(None, geneProdExcl_list))

    # Excluding some elements from gene/product name database    
    geneProdDB_list = list(set(geneProdDB_list)-set(keywVariants_list))
    geneProdDB_list = list(set(geneProdDB_list)-set(toxinDB_list))
    geneProdDB_list = list(set(geneProdDB_list)-set(geneProdExcl_list))
    
    output_fn = args.outfolder +\
                args.outfn_stem +\
                datetime.datetime.today().strftime('%Y_%m_%d_%H%M')

    ### STEP 2. Set up logger
    log = logging.getLogger(__name__)
    if verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.DEBUG, logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level=logging.INFO, logger=log)

    
    ### STEP 3. Setting up masterDict, querying PubMed
    masterDict = {}
    masterDict['pubmed_query'] = {}
    masterDict['pubmed_query']['query_string'] = re.sub('\s+',' ', query)
    ## Querying PubMed for user-supplied query
    action = "querying PubMed for user-supplied query"
    log.info("%s" % action)
    try:
        masterDict['pubmed_query']['query_return'] = PubMedInteract(email).search_via_query(query)
        masterDict['pubmed_query']['query_return']['IdList']
    except:
        raise Exception("Error when %s" % action)
    log.info("Number of uids matching the query found: %s" % len(masterDict['pubmed_query']['query_return']['IdList']))
    ## Updating output file
    action = "generating JSON output file"
    log.info("%s" % action)
    MiscOps().write_to_json(output_fn+'.json', masterDict)


    ### STEP 4. Retrieving metadata from PubMed, extracting metadata
    ## Looping through query uids
    action = "looping through query uids"
    log.info("%s" % action)
    for uid in masterDict['pubmed_query']['query_return']['IdList']:
        ## Setting up entry in masterDict
        masterDict[uid] = {}
        #masterDict[uid]['article_uid'] = uid
        ## Retrieving metadata from PubMed
        action = "retrieving metadata from PubMed"
        log.info("\tuid %s: %s" % (uid, action))
        try:
            xml_result = PubMedInteract(email).retrieve_metadata(uid)
            article = xml_result['PubmedArticle'][0]
        except:
            log.critical("\tuid %s: Error when %s" % (uid, action))
        ## Parsing the metadata, extracting relevant info
        action = "parsing the metadata, extracting relevant info"
        log.info("\tuid %s: %s" % (uid, action))
        try:
            masterDict[uid]['article_title'] = XMLparsing().article_title(article)
            masterDict[uid]['article_publdate'] = XMLparsing().article_publdate(article)
            masterDict[uid]['article_publdate_str'] = MiscOps().convert_to_date(masterDict[uid]['article_publdate'])
            masterDict[uid]['article_keywords'] = XMLparsing().article_keywords(article)
            masterDict[uid]['article_mesh'] = XMLparsing().article_mesh(article)
            masterDict[uid]['article_pmid'] = XMLparsing().article_pmid(article)
            masterDict[uid]['article_pmcid'] = PMCIDops().article_pmcid(masterDict[uid]['article_pmid'], email)
        except:
            log.critical("\tuid %s: Error when %s" % (uid, action))

    ## Updating output file
    action = "updating JSON output file"
    log.info("%s" % action)
    MiscOps().write_to_json(output_fn+'.json', masterDict)


    ### STEP 5. Extract and parse the fulltext of target publications
    action = "looping through masterDict items"
    log.info("%s" % action)
    for uid, article in masterDict.items():
        if uid == 'pubmed_query':
            continue
        if 'article_pmcid' in article.keys() and article['article_pmcid']:
            
            ## Forming BioC PMCID URL and retrieving full articles in XML format via BioC API
            action = "retrieving complete article from PubMedCentral"
            log.info("\tuid %s: %s" % (uid, action))
            try:
                article['article_bioc_url'] = PMCIDops().get_pmcid_bioc_url(article['article_pmcid'])
                article['article_pubreader_url'] = PMCIDops().get_pmcid_pubreader_url(article['article_pmcid'])
                try:
                    handle = urllib.request.urlopen(article['article_bioc_url'])
                    article_complete = handle.read()
                except:
                    try:
                        handle = urllib.request.Request(article['article_pubreader_url'], headers={'User-Agent': 'Mozilla/5.0'})
                        article_complete = urllib.request.urlopen(handle).read()
                    except:
                        log.warning("\tuid %s: Retrieval unsuccessful" % uid)
                        article_complete = None
                        article['article_paragraphs_with_keyw'] = None
            except:
                log.critical("\tuid %s: Error when %s" % (uid, action))

            if article_complete:
                ## Extracting the full text
                action = "extracting full text"
                log.info("\tuid %s: %s" % (uid, action))
                try:
                    if '<?xml version="1.0"' in str(article_complete):
                        # Parse the XML of the full text
                        fulltext_etree = lxml.etree.fromstring(article_complete)
                        # Remove unnecessary sections from full text
                        LXMLops(fulltext_etree).remove_expendable()
                        # Extract all text of the full text by paragraph
                        all_paragraphs = LXMLops(fulltext_etree).extract_all_text()
                        masterDict[uid]['n_fulltext_paragr'] = len(all_paragraphs)
                    if '<!DOCTYPE html>' in str(article_complete):
                        # Parse the HTML of the full text
                        fulltext_soup = bs4.BeautifulSoup(article_complete, 'html.parser')
                        # Remove unnecessary sections from full text
                        html_divs = HTMLops().remove_expendable(fulltext_soup)
                        # Extract all text of the full text by paragraph
                        all_paragraphs = HTMLops().extract_all_text(html_divs)
                        masterDict[uid]['n_fulltext_paragr'] = len(all_paragraphs)
                except:
                    log.critical("\tuid %s: Error when %s" % (uid, action))
                    
                ## PARAGRAPH FILTERING
                ## Extracting those paragraphs that contain the keyword `gene` in its different forms
                action = "extracting paragraphs with keyword `gene`"
                log.info("\tuid %s: %s" % (uid, action))
                try:
                    paragraphs_with_keyw = ParagrOps().extract_gene_paragraphs(all_paragraphs, keywVariants_list)
                    masterDict[uid]['article_paragraphs_with_keyw'] = paragraphs_with_keyw
                except:
                    log.critical("\tuid %s: Error when %s" % (uid, action))
        else:
            log.warning("\tuid %s: publication on PubMed but not on PubMedCentral" % uid)
            masterDict[uid]['article_paragraphs_with_keyw'] = None

    ## Updating output file
    action = "updating JSON output file"
    log.info("%s" % action)
    MiscOps().write_to_json(output_fn+'.json', masterDict)


    ### STEP 6. Extract and parse the fulltext of target publications
    ## Extract those paragraphs that match any of the the gene names of a gene name database
    action = "looping through masterDict items"
    log.info("%s" % action)
    for uid, article in masterDict.items():
        if uid == 'pubmed_query':
            continue
        if article['article_paragraphs_with_keyw']:
            
            ## PARAGRAPH SELECTION
            ## Identifying paragraphs that match any of the database names for both genes and toxins
            action = "identifying paragraphs that match any of the database names for both genes and toxins"
            log.info("\tuid %s: %s" % (uid, action))
            try:
                result_paragraphs, found_genes, found_toxins = [], [], []
                for paragr in article['article_paragraphs_with_keyw']:
                    for toxin in toxinDB_list:
                        if toxin.casefold() in paragr.casefold():
                            for gene in geneProdDB_list:
                                if gene.casefold() in paragr.casefold():
                                    result_paragraphs.append(paragr)
                                    found_toxins.append(toxin)
                                    found_genes.append(gene)
                if len(result_paragraphs) >= 1:
                    masterDict[uid]['article_result_paragraphs'] = ParagrOps().remove_duplicate_paragraphs(result_paragraphs)  # Note: The same paragraph may have been found multiple times.
                    masterDict[uid]['found_toxins'] = list(set(found_toxins))
                    masterDict[uid]['found_genes'] = list(set(found_genes))
                else:
                    log.info("\tuid %s: no fitting paragraphs detected" % uid)
                    masterDict[uid]['article_result_paragraphs'] = None
                    masterDict[uid]['found_toxins'] = None
                    masterDict[uid]['found_genes'] = None
            except:
                log.critical("\tuid %s: Error when %s" % (uid, action))
                    
        else:
            masterDict[uid]['article_result_paragraphs'] = None
            masterDict[uid]['found_toxins'] = None
            masterDict[uid]['found_genes'] = None

    ## Updating output file
    action = "updating JSON output file"
    log.info("%s" % action)
    MiscOps().write_to_json(output_fn+'.json', masterDict)


    ### STEP 7. Ordering articles by date of publication, with most recent publication first
    action = "ordering articles by date of publication, with most recent publication first"
    log.info("%s" % action)

    handleDict = {}
    tmpDict = dict(masterDict)  # Must copy masterDict; otherwise pop in tmpDict will also delete in masterDict
    handleDict['pubmed_query'] = tmpDict.pop('pubmed_query')
    handleDict.update(collections.OrderedDict(sorted(tmpDict.items(), key=lambda item: item[1]['article_publdate_str'], reverse=True)))
    masterDict = handleDict

    ## Updating output file
    action = "updating JSON output file"
    log.info("%s" % action)
    MiscOps().write_to_json(output_fn+'.json', masterDict)
    
#-----------------------------------------------------------------#
# MAIN

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Author|Version: '+__version__)#description="  --  ".join([__author__+' <'+__email__+'>', __info__, __version__]))
    parser.add_argument("--query", "-q", type=str, required=True, 
                        help="Query string to query NCBI PubMed via Entrez")
    parser.add_argument("--mail", "-m", type=str, required=True,
                        help="Your email address (needed for querying NCBI PubMed via Entrez)")
    parser.add_argument("--keywVariants_list", "-k", type=list, required=False, 
                        default=[' gene ', ' gene,', ' gene;', ' gene.'], 
                        help="(Optional) List of keyword variants")
    parser.add_argument("--toxinDB_fname", "-x", type=str, required=False, 
                        default="/home/mgruenst/tmp/databases/mycotoxin_names_2022_07_26.txt",
                        help="(Optional) Path to toxin name database file")
    parser.add_argument("--geneProdDB_fname", "-g", type=str, required=False, 
                        default="/home/mgruenst/tmp/databases/Trichoderma_genomes_GeneList_2022_07_27.txt",
                        help="(Optional) Path to gene/product database file")
    parser.add_argument("--geneProdExcl_fname", "-e", type=str, required=False, 
                        default="/home/mgruenst/tmp/databases/Trichoderma_genomes_GeneList_Exclusions_2022_07_27.txt",
                        help="(Optional) Path to file containing an list of gene/product names to be excluded")
    parser.add_argument("--outfn_stem", "-s", type=str, required=False, 
                        default="litMiningPubmed_Results_",
                        help="Stem of the output file")
    parser.add_argument("--outfolder", "-o", type=str, required=False, 
                        default="./results/",
                        help="Name of the output folder")
    parser.add_argument("--verbose", "-v", action="store_true", required=False, 
                        default=True, help="(Optional) Enable verbose logging")
    args = parser.parse_args()
    #if bool(args.query) ^ bool(args.mail):
    #    parser.error("--query and --mail must be given together")
    main(args)

#-----------------------------------------------------------------#
