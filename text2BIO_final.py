# -*- coding: utf-8 -*-
# @Author: Avijit Mitra
# Takes a text file, annotated xml file and create a dataset following BIO tagscheme


import sys
import argparse
import nltk
from nltk import word_tokenize, sent_tokenize
from nltk.tokenize import TreebankWordTokenizer as twt
import xml.etree.ElementTree as ET
import os.path


def convert(ehr_file, ann_file, output_file, limited, exclusion):
    '''Convert an EHR record to standard BIO scheme'''

    l_entities = ['BLEEDING_EVENT', 'DRUGNAME', 'TRIGGER_ALTERNATIVE_CAUSE', 'SEVERITY', 'BLEEDING_ANATOMIC_SITE', 'BLEEDING_LAB_EVAL']
    e_entities = ['NOT_USEFUL_FLAG']
    e_flag = False
    
    with open(ehr_file,'r') as ehr:
        ehr_text = ehr.read()
        # ehr_text = ehr_text.replace('=',' = ')
        # ehr_text = ehr_text.replace('-',' - ')
        # ehr_text = ehr_text.replace('/',' / ')
        # ehr_text = ehr_text.replace('~',' ~ ')
        # ehr_text = ehr_text.replace('+',' + ')
        tokenized_words, w_spans = [], []
        tokenizer = nltk.data.load('nltk:tokenizers/punkt/english.pickle')
        s_spans = [(s_ind, e_ind) for s_ind, e_ind in tokenizer.span_tokenize(ehr_text)]
            
        for index,s in enumerate(sent_tokenize(ehr_text)):
            tokenized_words += word_tokenize(s)
            w_spans += [(s_ind+s_spans[index][0],e_ind+s_spans[index][0]) for s_ind, e_ind in twt().span_tokenize(s)]
            tokenized_words += ['\n']                                           # all the words in the document, tokenized
            w_spans += ['\n']
            
    with open(ann_file,'r') as ann:
        ehr_tree = ET.parse(ann_file).getroot()
        labels, text_spans = [], []
        labels_dict, ind_dict = {}, {}
        for child in ehr_tree:
            if child.tag == 'classMention':
                for c in child:
                    if c.tag == 'mentionClass':
                        if limited and c.attrib['id'] not in l_entities: 
                            if c.attrib['id'] in e_entities: e_flag = True      # keep track if any e_entity is found
                            continue                                            # skip over the next ent. when an unimp. ent. is found 
                        labels += [c.attrib['id']]                              # entity as label
                        text_spans += [c.text]                                  # text span for the label
                        labels_dict[child.attrib['id']] = c.attrib['id']
            if child.tag == 'annotation':
                for c in child:
                    if c.tag == 'span':
                        s_ind = int(c.attrib['start'])
                        e_ind = int(c.attrib['end'])
                    if c.tag == 'mention':
                        entity_id = c.attrib['id']
                ind_dict[entity_id] = (s_ind,e_ind)
    
    if exclusion and e_flag:
        with open('excluded_files123.txt','a') as out:
            out.write(ehr_file+'\n')                                            # Store the file name , when ehr has e_entity
    else:
        with open(output_file,'a') as out:
            linecount = 0
            out.write('-DOCSTART- O\n\n')
            linecount += 2           
            for index,w in enumerate(tokenized_words):
                if w != '\n':
                    flag = 0
                    overlapped_b_entities, overlapped_i_entities, new_words, new_entities = [], [], [], []
                    for k,v in labels_dict.items():
                        if w_spans[index][0] ==  ind_dict[k][0]:                     # for the B- tags
                            overlapped_b_entities.append(v)
                            flag = 1
                            continue
                        elif ind_dict[k][0] < w_spans[index][0] < ind_dict[k][1]:    # for the I- tags
                            overlapped_i_entities.append(v)
                            flag = 1
                            continue
                        elif w_spans[index][0] < ind_dict[k][0] < w_spans[index][1]: # for annotation-tokenization mismatch, like '-HGB' and 'HGB'
                            flag = 2
                            new_words += [w[:ind_dict[k][0]-w_spans[index][0]], w[ind_dict[k][0]-w_spans[index][0]:]]
                            new_entities += ['O', 'B-'+v]
                        elif w_spans[index][0] < ind_dict[k][1] < w_spans[index][1]: # for annotation-tokenization mismatch, like '11.3*' and '11.3'
                            flag = 2
                            new_words += [w[:ind_dict[k][1]-w_spans[index][0]], w[ind_dict[k][1]-w_spans[index][0]:]]
                            if ind_dict[k][0] < w_spans[index][0]: new_entities += ['I-'+v, 'O']
                            else: new_entities += ['B-'+v, 'O']

                    # Resolve overlaps, I tags get priority over B tags.
                    if flag == 2:
                        for i,j in zip(new_words,new_entities): out.write(i+' '+j+'\n')
                    elif len(set(overlapped_i_entities)) == 1:
                        out.write(w + ' I-' + overlapped_i_entities[0] + '\n')
                        linecount += 1
                    elif 'BLEEDING_EVENT' in overlapped_i_entities:
                        out.write(w + ' I-BLEEDING_EVENT' + '\n')
                        linecount += 1
                    elif 'TRIGGER_ALTERNATIVE_CAUSE' in overlapped_i_entities:
                        out.write(w + ' I-TRIGGER_ALTERNATIVE_CAUSE' + '\n')
                        linecount += 1                        
                    elif 'BLEEDING_LAB_EVAL' in overlapped_i_entities and 'SEVERITY' in overlapped_i_entities:
                        out.write(w + ' I-BLEEDING_LAB_EVAL' + '\n')
                        linecount += 1
                    elif len(set(overlapped_b_entities)) == 1:
                        out.write(w + ' B-' + overlapped_b_entities[0] + '\n')
                        linecount += 1
                    elif 'BLEEDING_EVENT' in overlapped_b_entities:
                        out.write(w + ' B-BLEEDING_EVENT' + '\n')
                        linecount += 1
                    elif 'TRIGGER_ALTERNATIVE_CAUSE' in overlapped_b_entities:
                        out.write(w + ' B-TRIGGER_ALTERNATIVE_CAUSE' + '\n')
                        linecount += 1                                      
                    elif 'BLEEDING_LAB_EVAL' in overlapped_b_entities and 'SEVERITY' in overlapped_b_entities:
                        out.write(w + ' B-BLEEDING_LAB_EVAL' + '\n')
                        linecount += 1
                        
                    if flag == 0:
                        out.write(w + ' O\n')
                        linecount += 1
                else:
                    out.write('\n')
                    linecount += 1   
            
        # keep track of cumulative total lines and which files are they coming from
        # prev_lc = 0
        # if os.path.isfile('linecount_test_new.txt'):
        #      with open('linecount_test_new.txt','r') as lc:
        #         prev_lc = lc.readlines()[-1].split(',')[0]
        #         prev_lc = int(prev_lc)
        # with open('linecount_test_new.txt','a') as lc:
        #     lc.write(str(linecount+prev_lc)+', '+ann_file+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converting text to BIO tag scheme...')
    parser.add_argument('--input_ehr',  help='Input text file', default='None')
    parser.add_argument('--input_ann',  help='Input annotated xml file', default='None')
    parser.add_argument('--output', help='Output text file with BIO tag scheme', default='None')
    parser.add_argument('--limited', action='store_true', help='Uses a limited number of entities.')
    parser.add_argument('--exclusion', action='store_true',  help='Excludes EHRs with specific entities.')
    
    args = parser.parse_args()
    
    if args.input_ehr == 'None':
        print('Please define an EHR file as input and run the script again.')
        sys.exit()
    else: input_ehr = args.input_ehr
        
    if args.input_ann == 'None':
        print('Please define an annotated xml file and run the script again.')
        sys.exit()
    else: input_ann = args.input_ann
        
    if args.output == 'None':output_file = args.input_ehr+'.BIO.txt'
    else:output_file = args.output
    
    convert(input_ehr,input_ann,output_file,args.limited,args.exclusion)