# Complete and Wobble codon maps
# A codon map is a dict 

def codonmap(name, gencode):
  """Possible names: complete, wobble, stoprisk, synonymcomplete, synonymwobble"""
  if name == "complete":
    return Complete
  if name == "stoprisk":
    return StopRiskComplete
  if name == "wobble" and gencode == 1:
    return Wobble_1
  if name == "wobble" and gencode == 12:
    return Wobble_12
  if name == "wobble" and gencode == 26:
    return Wobble_26
  if name == "synonymcomplete" and gencode == 1:
    return SynonymComplete_1
  if name == "synonymcomplete" and gencode == 12:
    return SynonymComplete_12
  if name == "synonymcomplete" and gencode == 26:
    return SynonymComplete_26
  if name == "synonymwobble" and gencode == 1:
    return SynonymWobble_1
  if name == "synonymwobble" and gencode == 12:
    return SynonymWobble_12
  if name == "synonymwobble" and gencode == 26:
    return SynonymWobble_26
  if name == "synonymaa" and gencode == 1:
    return SynonymAA_1
  if name == "synonymaa" and gencode == 12:
    return SynonymAA_12
  if name == "synonymaa" and gencode == 26:
    return SynonymAA_26

# A Complete mapping is the identity mapping
Complete = {
'AAA':('AAA'), 'AAC':('AAC'), 'AAG':('AAG'), 'AAU':('AAU'),
'ACA':('ACA'), 'ACC':('ACC'), 'ACG':('ACG'), 'ACU':('ACU'),
'AGA':('AGA'), 'AGC':('AGC'), 'AGG':('AGG'), 'AGU':('AGU'),
'AUA':('AUA'), 'AUC':('AUC'), 'AUG':('AUG'), 'AUU':('AUU'),
'CAA':('CAA'), 'CAC':('CAC'), 'CAG':('CAG'), 'CAU':('CAU'),
'CCA':('CCA'), 'CCC':('CCC'), 'CCG':('CCG'), 'CCU':('CCU'),
'CGA':('CGA'), 'CGC':('CGC'), 'CGG':('CGG'), 'CGU':('CGU'),
'CUA':('CUA'), 'CUC':('CUC'), 'CUG':('CUG'), 'CUU':('CUU'),       
'GAA':('GAA'), 'GAC':('GAC'), 'GAG':('GAG'), 'GAU':('GAU'),
'GCA':('GCA'), 'GCC':('GCC'), 'GCG':('GCG'), 'GCU':('GCU'),
'GGA':('GGA'), 'GGC':('GGC'), 'GGG':('GGG'), 'GGU':('GGU'),
'GUA':('GUA'), 'GUC':('GUC'), 'GUG':('GUG'), 'GUU':('GUU'),
'UAA':('UAA'), 'UAC':('UAC'), 'UAG':('UAG'), 'UAU':('UAU'),
'UCA':('UCA'), 'UCC':('UCC'), 'UCG':('UCG'), 'UCU':('UCU'),
'UGA':('UGA'), 'UGC':('UGC'), 'UGG':('UGG'), 'UGU':('UGU'),
'UUA':('UUA'), 'UUC':('UUC'), 'UUG':('UUG'), 'UUU':('UUU')
}

# A Wobble mapping has codons recognized by the same tRNA
# use with genetic code 1
Wobble_1 = {
'AAA':('AAA'), 'AAC':('AAY'), 'AAG':('AAG'), 'AAU':('AAY'),
'ACA':('ACA'), 'ACC':('ACY'), 'ACG':('ACG'), 'ACU':('ACY'),
'AGA':('AGA'), 'AGC':('AGY'), 'AGG':('AGG'), 'AGU':('AGY'),
'AUA':('AUA'), 'AUC':('AUY'), 'AUG':('AUG'), 'AUU':('AUY'),
'CAA':('CAA'), 'CAC':('CAY'), 'CAG':('CAG'), 'CAU':('CAY'),
'CCA':('CCA'), 'CCC':('CCY'), 'CCG':('CCG'), 'CCU':('CCY'),
'CGA':('CGH'), 'CGC':('CGH'), 'CGG':('CGG'), 'CGU':('CGH'),
'CUA':('CUA'), 'CUC':('CUY'), 'CUG':('CUG'), 'CUU':('CUY'), 
'GAA':('GAA'), 'GAC':('GAY'), 'GAG':('GAG'), 'GAU':('GAY'),
'GCA':('GCA'), 'GCC':('GCY'), 'GCG':('GCG'), 'GCU':('GCY'),
'GGA':('GGA'), 'GGC':('GGY'), 'GGG':('GGG'), 'GGU':('GGY'),
'GUA':('GUA'), 'GUC':('GUY'), 'GUG':('GUG'), 'GUU':('GUY'),
'UAA':('UAA'), 'UAC':('UAY'), 'UAG':('UAG'), 'UAU':('UAY'),
'UCA':('UCA'), 'UCC':('UCY'), 'UCG':('UCG'), 'UCU':('UCY'),
'UGA':('UGA'), 'UGC':('UGY'), 'UGG':('UGG'), 'UGU':('UGY'),
'UUA':('UUA'), 'UUC':('UUY'), 'UUG':('UUG'), 'UUU':('UUY')
}

# Wobble mapping for genetic code 12
Wobble_12 = {
'AAA':('AAA'), 'AAC':('AAY'), 'AAG':('AAG'), 'AAU':('AAY'),
'ACA':('ACA'), 'ACC':('ACY'), 'ACG':('ACG'), 'ACU':('ACY'),
'AGA':('AGA'), 'AGC':('AGY'), 'AGG':('AGG'), 'AGU':('AGY'),
'AUA':('AUA'), 'AUC':('AUY'), 'AUG':('AUG'), 'AUU':('AUY'),
'CAA':('CAA'), 'CAC':('CAY'), 'CAG':('CAG'), 'CAU':('CAY'),
'CCA':('CCA'), 'CCC':('CCY'), 'CCG':('CCG'), 'CCU':('CCY'),
'CGA':('CGH'), 'CGC':('CGH'), 'CGG':('CGG'), 'CGU':('CGH'),
'CUA':('CUH'), 'CUC':('CUH'), 'CUG':('CUG'), 'CUU':('CUH'),
'GAA':('GAA'), 'GAC':('GAY'), 'GAG':('GAG'), 'GAU':('GAY'),
'GCA':('GCA'), 'GCC':('GCY'), 'GCG':('GCG'), 'GCU':('GCY'),
'GGA':('GGA'), 'GGC':('GGY'), 'GGG':('GGG'), 'GGU':('GGY'),
'GUA':('GUA'), 'GUC':('GUY'), 'GUG':('GUG'), 'GUU':('GUY'),
'UAA':('UAA'), 'UAC':('UAY'), 'UAG':('UAG'), 'UAU':('UAY'),
'UCA':('UCA'), 'UCC':('UCY'), 'UCG':('UCG'), 'UCU':('UCY'),
'UGA':('UGA'), 'UGC':('UGY'), 'UGG':('UGG'), 'UGU':('UGY'),
'UUA':('UUA'), 'UUC':('UUY'), 'UUG':('UUG'), 'UUU':('UUY')
}
  
# Wobble mapping for genetic code 12
Wobble_26 = {
'AAA':('AAA'), 'AAC':('AAY'), 'AAG':('AAG'), 'AAU':('AAY'),
'ACA':('ACA'), 'ACC':('ACY'), 'ACG':('ACG'), 'ACU':('ACY'),
'AGA':('AGA'), 'AGC':('AGY'), 'AGG':('AGG'), 'AGU':('AGY'),
'AUA':('AUA'), 'AUC':('AUY'), 'AUG':('AUG'), 'AUU':('AUY'),
'CAA':('CAA'), 'CAC':('CAY'), 'CAG':('CAG'), 'CAU':('CAY'),
'CCA':('CCA'), 'CCC':('CCY'), 'CCG':('CCG'), 'CCU':('CCY'),
'CGA':('CGH'), 'CGC':('CGH'), 'CGG':('CGG'), 'CGU':('CGH'),
'CUA':('CUH'), 'CUC':('CUH'), 'CUG':('CUG'), 'CUU':('CUH'),
'GAA':('GAA'), 'GAC':('GAY'), 'GAG':('GAG'), 'GAU':('GAY'),
'GCA':('GCA'), 'GCC':('GCY'), 'GCG':('GCG'), 'GCU':('GCY'),
'GGA':('GGA'), 'GGC':('GGY'), 'GGG':('GGG'), 'GGU':('GGY'),
'GUA':('GUA'), 'GUC':('GUY'), 'GUG':('GUG'), 'GUU':('GUY'),
'UAA':('UAA'), 'UAC':('UAY'), 'UAG':('UAG'), 'UAU':('UAY'),
'UCA':('UCA'), 'UCC':('UCY'), 'UCG':('UCG'), 'UCU':('UCY'),
'UGA':('UGA'), 'UGC':('UGY'), 'UGG':('UGG'), 'UGU':('UGY'),
'UUA':('UUA'), 'UUC':('UUY'), 'UUG':('UUG'), 'UUU':('UUY')
}

# A complete stop risk map indicates how many single mutations
# can transform a codon into a stop codon, the codons are 
# defined as in Complete CodonMap
StopRiskComplete = {
'AAA':(1), 'AAC':(0), 'AAG':(1), 'AAU':(0),
'ACA':(0), 'ACC':(0), 'ACG':(0), 'ACU':(0),
'AGA':(1), 'AGC':(0), 'AGG':(0), 'AGU':(0),
'AUA':(0), 'AUC':(0), 'AUG':(0), 'AUU':(0),
'CAA':(1), 'CAC':(0), 'CAG':(1), 'CAU':(0),
'CCA':(0), 'CCC':(0), 'CCG':(0), 'CCU':(0),
'CGA':(1), 'CGC':(0), 'CGG':(0), 'CGU':(0),
'CUA':(0), 'CUC':(0), 'CUG':(0), 'CUU':(0),
'GAA':(1), 'GAC':(0), 'GAG':(1), 'GAU':(0),
'GCA':(0), 'GCC':(0), 'GCG':(0), 'GCU':(0),
'GGA':(1), 'GGC':(0), 'GGG':(0), 'GGU':(0),
'GUA':(0), 'GUC':(0), 'GUG':(0), 'GUU':(0),
'UAA':(0), 'UAC':(2), 'UAG':(0), 'UAU':(2),
'UCA':(2), 'UCC':(0), 'UCG':(1), 'UCU':(0),
'UGA':(0), 'UGC':(1), 'UGG':(2), 'UGU':(1),
'UUA':(2), 'UUC':(0), 'UUG':(1), 'UUU':(0)
}

# A complete synonym map indicates the list of codons which
# are translated into the same amino-acid, with the same 
# first and second nucleotides. There are several 
# synonym maps according to the genetic code used for 
# translation

SynonymComplete_1 = {
'AAR':('AAA','AAG'),
'AAY':('AAC','AAU'),
'ACN':('ACA','ACC','ACG','ACU'),
'AGR':('AGA','AGG'),
'AGY':('AGC','AGU'),
'AUG':('AUG',),
'AUH':('AUA','AUC','AUU'),
'CAR':('CAA','CAG'),
'CAY':('CAC','CAU'),
'CCN':('CCA','CCC','CCG','CCU'),
'CGN':('CGA','CGC','CGG','CGU'),
'CUN':('CUA','CUC','CUG','CUU'),
'GAR':('GAA','GAG'),
'GAY':('GAC','GAU'),
'GCN':('GCA','GCC','GCG','GCU'),
'GGN':('GGA','GGC','GGG','GGU'),
'GUN':('GUA','GUC','GUG','GUU'),
'UAY':('UAC','UAU'),
'UCN':('UCA','UCC','UCG','UCU'),
'UGG':('UGG',),
'UGY':('UGC','UGU'),
'UUR':('UUA','UUG'),
'UUY':('UUC','UUU'),
}

SynonymComplete_12 = {
'AAR':('AAA','AAG'),
'AAY':('AAC','AAU'),
'ACN':('ACA','ACC','ACG','ACU'),
'AGR':('AGA','AGG'),
'AGY':('AGC','AGU'),
'AUG':('AUG',),
'AUH':('AUA','AUC','AUU'),
'CAR':('CAA','CAG'),
'CAY':('CAC','CAU'),
'CCN':('CCA','CCC','CCG','CCU'),
'CGN':('CGA','CGC','CGG','CGU'),
'CUG':('CUG',),
'CUH':('CUA','CUC','CUU'),
'GAR':('GAA','GAG'),
'GAY':('GAC','GAU'),
'GCN':('GCA','GCC','GCG','GCU'),
'GGN':('GGA','GGC','GGG','GGU'),
'GUN':('GUA','GUC','GUG','GUU'),
'UAY':('UAC','UAU'),
'UCN':('UCA','UCC','UCG','UCU'),
'UGG':('UGG',),
'UGY':('UGC','UGU'),
'UUR':('UUA','UUG'),
'UUY':('UUC','UUU'),
}

SynonymComplete_26 = {
'AAR':('AAA','AAG'),
'AAY':('AAC','AAU'),
'ACN':('ACA','ACC','ACG','ACU'),
'AGR':('AGA','AGG'),
'AGY':('AGC','AGU'),
'AUG':('AUG',),
'AUH':('AUA','AUC','AUU'),
'CAR':('CAA','CAG'),
'CAY':('CAC','CAU'),
'CCN':('CCA','CCC','CCG','CCU'),
'CGN':('CGA','CGC','CGG','CGU'),
'CUG':('CUG',),
'CUH':('CUA','CUC','CUU'),
'GAR':('GAA','GAG'),
'GAY':('GAC','GAU'),
'GCN':('GCA','GCC','GCG','GCU'),
'GGN':('GGA','GGC','GGG','GGU'),
'GUN':('GUA','GUC','GUG','GUU'),
'UAY':('UAC','UAU'),
'UCN':('UCA','UCC','UCG','UCU'),
'UGG':('UGG',),
'UGY':('UGC','UGU'),
'UUR':('UUA','UUG'),
'UUY':('UUC','UUU'),
}

SynonymWobble_1 = {
'AAR':('AAA','AAG'),
'AAY':('AAY',),
'ACN':('ACA','ACG','ACY'),
'AGR':('AGA','AGG'),
'AGY':('AGY',),
'AUG':('AUG',),
'AUH':('AUA','AUY'),
'CAR':('CAA','CAG'),
'CAY':('CAY',),
'CCN':('CCA','CCG','CCY'),
'CGN':('CGH','CGG'),
'CUN':('CUA','CUG','CUY'),
'GAR':('GAA','GAG'),
'GAY':('GAY',),
'GCN':('GCA','GCG','GCY'),
'GGN':('GGA','GGG','GGY'),
'GUN':('GUA','GUG','GUY'),    
'UAY':('UAY',),
'UCN':('UCA','UCG','UCY'),
'UGG':('UGG',),
'UGY':('UGY',),
'UUR':('UUA','UUG'),
'UUY':('UUY',),
}

SynonymWobble_12 = {
'AAR':('AAA','AAG'),
'AAY':('AAY',),
'ACN':('ACA','ACG','ACY'),
'AGR':('AGA','AGG'),
'AGY':('AGY',),
'AUG':('AUG',),
'AUH':('AUA','AUY'),
'CAR':('CAA','CAG'),
'CAY':('CAY',),
'CCN':('CCA','CCG','CCY'),
'CGN':('CGH','CGG'),
'CUG':('CUG',),
'CUH':('CUH',),
'GAR':('GAA','GAG'),
'GAY':('GAY',),
'GCN':('GCA','GCG','GCY'),
'GGN':('GGA','GGG','GGY'),
'GUN':('GUA','GUG','GUY'),
'UAY':('UAY',),
'UCN':('UCA','UCG','UCY'),
'UGG':('UGG',),
'UGY':('UGY',),
'UUR':('UUA','UUG'),
'UUY':('UUY',),
}

SynonymWobble_26 = {
'AAR':('AAA','AAG'),
'AAY':('AAY',),
'ACN':('ACA','ACG','ACY'),
'AGR':('AGA','AGG'),
'AGY':('AGY',),
'AUG':('AUG',),
'AUH':('AUA','AUY'),
'CAR':('CAA','CAG'),
'CAY':('CAY',),
'CCN':('CCA','CCG','CCY'),
'CGN':('CGH','CGG'),
'CUG':('CUG',),
'CUH':('CUH',),
'GAR':('GAA','GAG'),
'GAY':('GAY',),
'GCN':('GCA','GCG','GCY'),
'GGN':('GGA','GGG','GGY'),
'GUN':('GUA','GUG','GUY'),
'UAY':('UAY',),
'UCN':('UCA','UCG','UCY'),
'UGG':('UGG',),
'UGY':('UGY',),
'UUR':('UUA','UUG'),
'UUY':('UUY',),
}

SynonymAA_1 = {
'Lys':('AAA','AAG'),
'Asn':('AAC','AAU'),
'Thr':('ACA','ACC','ACG','ACU'),
'Met':('AUG',),
'Ile':('AUA','AUC','AUU'),
'Gln':('CAA','CAG'),
'His':('CAC','CAU'),
'Pro':('CCA','CCC','CCG','CCU'),
'Arg':('CGA','CGC','CGG','CGU','AGA','AGG'),
'Glu':('GAA','GAG'),
'Asp':('GAC','GAU'),
'Ala':('GCA','GCC','GCG','GCU'),
'Gly':('GGA','GGC','GGG','GGU'),
'Val':('GUA','GUC','GUG','GUU'),
'Tyr':('UAC','UAU'),
'Ser':('UCA','UCC','UCG','UCU','AGC','AGU'),
'Trp':('UGG',),
'Cys':('UGC','UGU'),
'Leu':('UUA','UUG','CUA','CUC','CUG','CUU'),
'Phe':('UUC','UUU'),
}

SynonymAA_12 = {
'Lys':('AAA','AAG'),
'Asn':('AAC','AAU'),
'Thr':('ACA','ACC','ACG','ACU'),
'Met':('AUG',),
'Ile':('AUA','AUC','AUU'),
'Gln':('CAA','CAG'),
'His':('CAC','CAU'),
'Pro':('CCA','CCC','CCG','CCU'),
'Arg':('CGA','CGC','CGG','CGU','AGA','AGG'),
'Glu':('GAA','GAG'),
'Asp':('GAC','GAU'),
'Ala':('GCA','GCC','GCG','GCU'),
'Gly':('GGA','GGC','GGG','GGU'),
'Val':('GUA','GUC','GUG','GUU'),
'Tyr':('UAC','UAU'),
'Ser':('UCA','UCC','UCG','UCU','AGC','AGU','CUG'),
'Trp':('UGG',),
'Cys':('UGC','UGU'),
'Leu':('UUA','UUG','CUA','CUC','CUU'),
'Phe':('UUC','UUU'),
}

SynonymAA_26 = {
'Lys':('AAA','AAG'),
'Asn':('AAC','AAU'),
'Thr':('ACA','ACC','ACG','ACU'),
'Met':('AUG',),
'Ile':('AUA','AUC','AUU'),
'Gln':('CAA','CAG'),
'His':('CAC','CAU'),
'Pro':('CCA','CCC','CCG','CCU'),
'Arg':('CGA','CGC','CGG','CGU','AGA','AGG'),
'Glu':('GAA','GAG'),
'Asp':('GAC','GAU'),
'Ala':('GCA','GCC','GCG','GCU','CUG'),
'Gly':('GGA','GGC','GGG','GGU'),
'Val':('GUA','GUC','GUG','GUU'),
'Tyr':('UAC','UAU'),
'Ser':('UCA','UCC','UCG','UCU','AGC','AGU'),
'Trp':('UGG',),
'Cys':('UGC','UGU'),
'Leu':('UUA','UUG','CUA','CUC','CUU'),
'Phe':('UUC','UUU'),
}

