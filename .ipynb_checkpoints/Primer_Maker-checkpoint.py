
import argparse

import pandas as pd
import numpy as np

from itertools import product

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import re

parser=argparse.ArgumentParser(description="Create primers from Sequences")
parser.add_argument('-i', "--inp", metavar="csv_file", help="Sequence input file" , type=str, required=True)
parser.add_argument('-o', "--out", metavar="filename", help="Primer List output filename" , type=str, required=True)
parser.add_argument("--well96", help="add this option to output in 96 well plate format (ensure you have <= 48 PCRs)", action='store_true')
parser.add_argument("--name96", help="name of 96 well plate", default="my_plate")
parser.add_argument('--gg5', metavar="(ATGC)x4", help="5 prime golden gate overhang. GG function assumes PCR is for scar-less insert")
parser.add_argument('--gg3', metavar="(ATGC)x4", help="3 prime golden gate 'underhang' (type out 5' to 3')")
args=parser.parse_args()


# Tm fxn    
def mytm(seq):
	return mt.Tm_NN(seq, dnac1=500, dnac2=0, selfcomp=False, Na=0, K=50, Tris=25, Mg=2, dNTPs=0.2, saltcorr=5)

# Fxn Tms for fwd sequence   
def tmfwd(dna):
    fwddna_list = []
    fwdtms_list = []
    count = 14
    while count < 61:
        count += 1
        if mytm(dna[0:count]) > 55 and mytm(dna[0:count]) < 72:
                fwddna_list.append(dna[0:count])
                fwdtms_list.append(mytm(dna[0:count]))
    seqcolm = pd.Series(fwddna_list,name='Sequence')
    tmcolm = pd.Series(fwdtms_list,name='Tms').apply(np.round).astype(int)
    return pd.concat([seqcolm,tmcolm],axis=1)

# Fxn Tms for rev 
def tmrev(dna): 
    seqdnas = Seq(dna)
    revcompl = str(seqdnas.reverse_complement())
    revdna_list = []
    revtms_list = []
    count = 14
    while count < 61:
        count += 1
        if mytm(dna[0:count]) > 55 and mytm(dna[0:count]) < 72:
                revdna_list.append(dna[0:count])
                revtms_list.append(mytm(dna[0:count]))
    revseqcolm = pd.Series(revdna_list,name='Sequence')
    revtmcolm = pd.Series(revtms_list,name='Tms').apply(np.round).astype(int)
    return pd.concat([revseqcolm,revtmcolm],axis=1)

# Main fxn
def primerfxn(dna):
    df3 = tmrev(dna) 
    Tmsrev = df3['Tms'].to_list()
    df2 = tmfwd(dna)
    Tmsfwd = df2['Tms'].to_list()
    # more definitions
    allposs = list(product(Tmsfwd,Tmsrev))
    finddiff = []
    for pair in allposs:
        finddiff.append(abs(pair[0]-pair[1]))
    findmindiff = allposs[np.argmin(finddiff)]
    fwdprim_list = (df2[df2['Tms'] == findmindiff[0]]['Sequence'].to_list()[0])
    fwdprimstms_list = (df2[df2['Tms'] == findmindiff[0]]['Tms'].to_list()[0])
    revprim_list = (df3[df3['Tms'] == findmindiff[1]]['Sequence'].to_list()[0])
    revprimstms_list = (df3[df3['Tms'] == findmindiff[1]]['Tms'].to_list()[0])
    return fwdprim_list, revprim_list ,revprimstms_list, fwdprimstms_list



#Anneal temp fxn
def get_anneal(temp1, temp2):
	return min(temp1, temp2)+1


# Primer maker for List of DNA 
def PrimerMaker(mydnalist):
    results = []
    for dnas in mydnalist: 
        results.append(primerfxn(dnas))
        df = pd.DataFrame(results, columns = ['fwd_primer' , 'rev_primer', 'Fwd Primer Tm', 'Rev Primer Tm'])
        df["anneal"] = df.apply(lambda x: get_anneal(x["Fwd Primer Tm"], x["Rev Primer Tm"]), axis=1)
    return df

def label_maker(gene_name, anneal, fwd=True, custom_str=""):
    if fwd == True:
        return gene_name + "_F_" + str(anneal) + custom_str
    if fwd == False:
        return gene_name + "_R_" + str(anneal) + custom_str

def series_joiner(ser1, ser2, ser3):
    return pd.concat([ser1, ser2, ser3], axis=1)

def bsai_regex(seq):
    my_matches = []
    match = re.finditer(r"(GGTCTC)|(GAGACC)", seq, flags=re.IGNORECASE) 
    for m in match:
        my_matches.append(m.span())
    return my_matches
def bsmbi_regex(seq):
    my_matches = []
    match = re.finditer(r"(CGTCTC)|(GAGACG)", seq, flags=re.IGNORECASE) 
    for m in match:
        my_matches.append(m.span())
    return my_matches 
def bbsi_regex(seq):
    my_matches = []
    match = re.finditer(r"(GAAGAC)|(GTCTTC)", seq, flags=re.IGNORECASE) 
    for m in match:
        my_matches.append(m.span())
    return my_matches
def bad_base(seq):
    my_matches = []
    match = re.finditer(r"[^ATGC]", seq, flags=re.IGNORECASE) 
    for m in match:
        my_matches.append(m.span())
    return my_matches

def recommender(l_bsai, l_bsmbi, l_bbsi):
    if len(l_bsai) == 0:
        return "BsaI"
    if len(l_bbsi) == 0:
        return "BbsI"
    if len(l_bsmbi) == 0:
        return "BsmbI"
    else:
        return "not good"
    
def type2s_wrapper(seq):
    bsai_list = bsai_regex(seq)
    bsmbi_list = bsmbi_regex(seq)
    bbsi_list = bbsi_regex(seq)
    return recommender(bsai_list, bsmbi_list, bbsi_list)

BsaIer = lambda x, y: "GGTCTCg" + y + x 
BsmbIer = lambda x, y: "CGTCTCg" + y + x 
BbsIer = lambda x, y: "GAAGACtg" + y + x 

def type2s_appender(seq, recommend, overhang):
    if recommend == "BsaI":
        return BsaIer(seq, overhang)
    if recommend == "BsmbI":
        return BsmbIer(seq, overhang)
    if recommend == "BbsI":
        return BbsIer(seq, overhang)

def char_range(c1, c2):
    """Generates the characters from `c1` to `c2`, inclusive."""
    for c in range(ord(c1), ord(c2)+1):
        yield chr(c)
    
def generate_96_labs():
    for number in range(1,13):
        for letter in char_range("A", "H"):
            yield(letter+str(number))
         
        
def append_well_labs(fwd_frame, rev_frame):
    fwd_frame.rename(columns={"label": "Name", "seq": "Sequence"}, inplace=True)
    rev_frame.rename(columns={"label": "Name", "seq": "Sequence"}, inplace=True)
    my_lab = generate_96_labs()
    for idx, entry in fwd_frame.iterrows():
        curr_lab = fwd_frame.at[idx, "Well Position"] = next(my_lab)
    while curr_lab[0] != "H":
        curr_lab = next(my_lab)
    for idx, entry in rev_frame.iterrows():
        rev_frame.at[idx, "Well Position"] = next(my_lab)
    result = pd.concat([fwd_frame, rev_frame])
    return result[result.columns[::-1]]

if args.gg5 is not None or args.gg3 is not None:
    assert len(args.gg5) == 4, "gg5 should be 4 bases"
    assert len(args.gg3) == 4, "gg3 should be 4 bases"


gene_df = pd.read_csv(args.inp)
gene_df.rename(columns={gene_df.columns[0]: "Name"}, inplace = True)
names = gene_df.iloc[:,0]
sequences = gene_df.iloc[:,1].to_list()
primer_df = PrimerMaker(sequences)

if args.gg5 is None or args.gg3 is None:
    result = pd.concat([names, primer_df[["fwd_primer", "rev_primer", "anneal"]]], axis=1)
    result["fwd_label"] = result.apply(lambda x: label_maker(x["Name"], x["anneal"]), axis=1)
    result["rev_label"] = result.apply(lambda x: label_maker(x["Name"], x["anneal"], fwd=False), axis=1)

    fwd_order = pd.DataFrame()
    fwd_order["label"] = result["fwd_label"] 
    fwd_order["seq"] = result["fwd_primer"]

    rev_order = pd.DataFrame()
    rev_order["label"] = result["rev_label"] 
    rev_order["seq"] = result["rev_primer"]
else:
    gene_df["recommend"] = gene_df.iloc[:,1].apply(type2s_wrapper)
    result = pd.concat([names, gene_df["recommend"], primer_df[["fwd_primer", "rev_primer", "anneal"]]], axis=1)
    result["fwd_label"] = result.apply(lambda x: label_maker(x["Name"], x["anneal"]), axis=1)
    result["rev_label"] = result.apply(lambda x: label_maker(x["Name"], x["anneal"], fwd=False), axis=1)

    fwd_order = pd.DataFrame()
    fwd_order["label"] = result["fwd_label"] 
    fwd_order["seq"] = result["fwd_primer"]
    fwd_order["recommend"] = result["recommend"]

    rev_order = pd.DataFrame()
    rev_order["label"] = result["rev_label"] 
    rev_order["seq"] = result["rev_primer"]
    rev_order["recommend"] = result["recommend"]

    fwd_order.seq = fwd_order.apply(lambda x: type2s_appender(x["seq"], x["recommend"], args.gg5), axis=1)
    rev_order.seq = rev_order.apply(lambda x: type2s_appender(x["seq"], x["recommend"], args.gg3), axis=1)

    fwd_order.dropna(inplace=True)
    fwd_order.drop(columns="recommend", inplace=True)
    fwd_order.reset_index(inplace=True, drop=True)

    rev_order.dropna(inplace=True)
    rev_order.drop(columns="recommend", inplace=True)
    rev_order.reset_index(inplace=True, drop=True)

pretty_print = pd.concat([fwd_order,rev_order])
pretty_print.reset_index(inplace=True, drop=True)
    
if args.well96:
    idt_order = append_well_labs(fwd_order, rev_order)
    idt_order.reset_index(inplace=True, drop=True)
else:    
    idt_order = pd.concat([fwd_order,rev_order])
    idt_order.reset_index(inplace=True, drop=True)

base_check = idt_order.iloc[:,1].apply(bad_base)
base_check_count = idt_order.iloc[:,1].apply(bad_base).apply(len).to_list()
bad_base_list = []
for idx, count in enumerate(base_check_count):
    if count != 0:
        bad_base_list.append(idx+1)
if len(bad_base_list) == 0:
    print("***Copy-paste into bulk input or open the excel file***\n")
    for idx, item in pretty_print.iterrows():
        print(item[0]+","+item[1])
else:
    print("Error! Ambiguous bases found at seq numbers: ", bad_base_list)
assert bad_base_list != 0, "Ambiguous bases found, please correct"
idt_order.to_excel(args.out+".xls", index=False, sheet_name=args.name96)