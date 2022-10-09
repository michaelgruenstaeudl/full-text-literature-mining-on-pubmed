Full text literature mining on PubMed for keyword combinations
==============================================================

A Python package for literature mining on PubMed

## EXAMPLE USAGE 1 - Toxins of Cyanobacteria

#### Cyanobacteria - conservative search[^1]
```
python3 01_litMiningPubmed_conductMining.py \
-q "(Cyanobacteria[tiab] OR Cyanobacteria[mh]) \
AND (cyanotox*[tiab] OR cyanotox*[mh]) \
AND (*synthesis[tiab] OR *synthesis[mh]) \
AND (gene[tiab] OR gene[mh]) \
AND 2000:2022[dp] \
AND free full text[sb]" \
-m your_email@address_here.at
```

#### Cyanobacteria - liberal search -- NOT RECOMMENDED
```
python3 01_litMiningPubmed_conductMining.py \
-q "Cyanobacteria[all]\
AND cyanotox*[all] \
AND *synthesis[all] \
AND gene[all] \
AND 2000:2022[dp] \
AND free full text[sb]" \
-m your_email@address_here.at
```


[^1]: Notes regarding query string: 
- "cyanotox*" is necessary b/c possible terms include "cyanotoxin(s)", "cyanotoxicity", etc.
- "2000:2022[dp]" selects only those publications that have been published in the genomic era (2000 ff.)
- "free full text[sb]" is necessary to ensure that full text of article can be retrieved
- "AND open access[filter]" would find all PMC Open Access subset articles with this filter but is too restrictive in practice


## EXAMPLE USAGE 2 - Mycotoxins of Trichoderma

#### Trichoderma - conservative search[^1]
```
python3 01_litMiningPubmed_conductMining.py \
-q "(Trichoderma[tiab] OR Trichoderma[mh]) \
AND (mycotox*[tiab] OR mycotox*[mh]) \
AND (*synthesis[tiab] OR *synthesis[mh]) \
AND (gene[tiab] OR gene[mh]) \
AND 2000:2022[dp] \
AND free full text[sb]" \
-m your_email@address_here.at
```

