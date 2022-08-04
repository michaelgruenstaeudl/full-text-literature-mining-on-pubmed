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
