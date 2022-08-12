Full text literature mining on PubMed for keyword combinations
==============================================================

A Python package for literature mining on PubMed

## EXAMPLE USAGE

#### Cyanobacteria - conservative search [^1]
```
python3 01_litMiningPubmed_conductMining.py \
-q "(Cyanobacteria[tiab] OR Cyanobacteria[mh]) \
AND (cyanotox*[tiab] OR cyanotox*[mh]) \
AND (*synthesis[tiab] OR *synthesis[mh]) \
AND (gene[tiab] OR gene[mh]) \
AND 2000:2022[dp] \
AND free full text[sb]" \
-m info@michael-gruenstaeudl.at
```

[^1]: Notes regarding query string: 
- "cyanotox*" is necessary b/c possible terms include "cyanotoxin(s)", "cyanotoxicity", etc.
- "2000:2022[dp]" selects only those publications that have been published in the genomic era (2000 ff.)
- "free full text[sb]" is necessary to ensure that full text of article can be retrieved
- "AND open access[filter]" would find all PMC Open Access subset articles with this filter but is too restrictive in practice
