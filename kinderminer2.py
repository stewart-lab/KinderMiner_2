import sys
import json
import datetime
import requests
import argparse
import os.path
import time

import lbd_stew as lbd

from dateutil.relativedelta import *

def build_arg_parser():
    parser = argparse.ArgumentParser(description='Run CHTC query')
    parser.add_argument('-s', '--sep', action='store_true', 
            help='perform keyphrase match as separate tokens')
    parser.add_argument('-t', '--targ_sep', action='store_true', 
            help='perform target term match as separate tokens')
    parser.add_argument('-e', '--eng', action='store_true', 
            help='perform keyphrase match with english language stemming')
    parser.add_argument('-f', '--targ_eng', action='store_true', 
            help='perform target term match with english language stemming')
    parser.add_argument('-a', '--alias', action='store_true', default=True,
            help='perform alias matching for target terms and keyphrase' +
            '(aliases for both use a specified delimiter [default:pipeline])')  
    parser.add_argument('-d', '--delimiter', default='|', 
            help='defines the delimiter to be used for alias separation' + 
            'when that option is used [default: pipeline]')
    parser.add_argument('-y', '--year', type=int, 
            default=datetime.datetime.now().year, 
            help='limit search to publications through particular year')
    parser.add_argument('term_file', 
            help='file containing all target terms to rank, one per line')
    parser.add_argument('keyphrase', 
            help='key phrase to rank target terms against')
    parser.add_argument('-db', '--db_version', type=str, default = None,
            help='PubMed Database version')
    parser.add_argument('-kf', '--keyphrasefile', action='store_true') 
    parser.add_argument('-ti', '--title_only', default=False, 
            action='store_true', help='query only in title')
    parser.add_argument('-ab', '--abstract_wo_abbr_expansion', 
            action='store_true',
            help='points to original PubMed database, with the abbrevation' +
            ' expansions turned off.' )
    parser.add_argument('-o', '--output_directory', type=str,
            default='',
            help='output directory')
    return parser

def check_args(args):
    if not args.keyphrase:
        print("Input Error: keyphrase file required")
        return False
    return True

def check_common_synonyms(key_phrase, target):
    kp_tokens_lowercase = [el.lower() for el in key_phrase.text]
    tt_tokens_lowercase = [el.lower() for el in target.text]
    if len(set(kp_tokens_lowercase).intersection(set(tt_tokens_lowercase))) > 0:
        return True
    return False

def get_output_file(key_phrase, out_dir):
    CUI_key = key_phrase.id
    kp_synonym_list = key_phrase.text
    cuikey = CUI_key.split('_')[0]
    syn = kp_synonym_list[0].replace(' ', '_')
    syn = syn.replace('/', '_')
    outfile_name = cuikey + '_' + syn
    outfile_path_name = os.path.join(out_dir, outfile_name + ".txt")
    return outfile_path_name

def check_output_file(outfile_path_name):
    all_lines = []
    targets_executed = []
    if os.path.isfile(outfile_path_name):
        for line in open(outfile_path_name):
            if line.startswith('target\t'):
                continue
            line_elements = line.strip().split('\t')
            if len(line_elements) == 5:
                all_lines.append(line.strip())
                targets_executed.append(line.strip().split(':')[0])
    return all_lines, targets_executed

def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    if not check_args(args):
        exit()

    # command line args
    TARGET_TERM_FILE = args.term_file
    KEY_PHRASE_FILE_NAME = args.keyphrase
    THROUGH_YEAR = args.year
    PUBMED_DB = args.db_version
    SEP_KP = args.sep
    STEM_KP = args.eng
    SEP_TARG = args.targ_sep
    STEM_TARG = args.targ_eng
    ALIAS = args.alias
    DELIM = args.delimiter
    TITLE = args.title_only
    ABSTRACT_IS_LONGFORM = not args.abstract_wo_abbr_expansion
    OUTPUT_DIR = args.output_directory

    # URL_BASE with user defined PubMed version
    URL_BASE = lbd.get_url_base(PUBMED_DB)

    # OUTPUT_DIR
    if OUTPUT_DIR:
        out_dir = OUTPUT_DIR

    # call chtc_query
    with open(KEY_PHRASE_FILE_NAME) as infile:
        keyphrase_list = [lbd.MatchText(k, SEP_KP, STEM_KP, ALIAS, TITLE, 
                            DELIM, ABSTRACT_IS_LONGFORM) for k in infile]

    with open(TARGET_TERM_FILE) as infile:
        target_list = [lbd.MatchText(l, SEP_TARG, STEM_TARG, ALIAS, TITLE, 
                            DELIM, ABSTRACT_IS_LONGFORM) for l in infile] 
    
    outfile_path_name, output_array, id_synonyms = perform_chtc_query(
                                                                keyphrase_list, 
                                                                target_list,
                                                                THROUGH_YEAR, 
                                                                URL_BASE, 
                                                                out_dir)

    # write output to a file specific to keyphrase
    write_to_file(outfile_path_name, output_array)    

def perform_chtc_query(keyphrase_list, target_list, THROUGH_YEAR, 
                                                            URL_BASE, out_dir):
    id_synonyms = {}
    
    # compute the total number of articles in the database
    db_article_cnt = lbd.get_count(None, None, THROUGH_YEAR, URL_BASE)

    for key_phrase in keyphrase_list:
        output_array = []

        # get contents output file (if exists)
        outfile_path_name = get_output_file(key_phrase, out_dir)
        output_array, targets_executed = check_output_file(outfile_path_name)   
        
        # check whether the file contains same number of target output
        if len(target_list) == len(output_array):
            continue    
    
        # check whether the file contains all target output
        if len(set(target_list).intersection(set(output_array))) == 0:
            continue
                    
        # individual keyphrase count
        kp_cnt = lbd.get_count(None, key_phrase, THROUGH_YEAR, URL_BASE)
        
        for target in target_list:
            targ_with_kp_cnt = 0

            # if target is in incomplete output file, skip execution
            if target.id in targets_executed:
                continue;

            # if key_phrase and target_term are same, skip execution
            if check_common_synonyms(key_phrase, target):
                continue

            # first the individual target count
            targ_cnt = lbd.get_count(target, None, THROUGH_YEAR, URL_BASE)

            # combined count, only if target and key phrase are present
            if targ_cnt > 0 and kp_cnt > 0:
                # now do both key phrase and target
                targ_with_kp_cnt = lbd.get_count(target, key_phrase, 
                                                        THROUGH_YEAR, URL_BASE)

            key = target.id
            value = target.text
            id_synonyms[key] = value
            outstr = '{0}\t{1}\t{2}\t{3}\t{4}'.format(key+':'+value[0], +
                                                        targ_with_kp_cnt, +
                                                        targ_cnt, +
                                                        kp_cnt, +
                                                        db_article_cnt)
            # print(outstr)
            output_array.append(outstr)
        
        # return output for every key phrase
        return outfile_path_name, output_array, id_synonyms
           
def write_to_file(outfile_path_name, output_array):
    with open(outfile_path_name, 'w') as out_fh:
        out_fh.write('target\ttarget_with_keyphrase_count\t' +
                        'target_count\tkeyphrase_count\tdb_article_count\n')
        for each_line in output_array:
            out_fh.write(each_line + '\n')


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()

    #print(end - start)
