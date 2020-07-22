import sys
import datetime
import time
import scipy.stats as stats
import numpy as np
import math

class KinderMinerCount(object):
    def __init__(self, target, target_and_keyphrase_count, target_count, 
                                        keyphrase_count, db_article_count):
        self.target = target
        self.target_and_keyphrase_count = target_and_keyphrase_count
        self.target_count = target_count
        self.keyphrase_count = keyphrase_count
        self.db_article_count = db_article_count
        self.notarget_count = db_article_count - target_count
        self.nokeyphrase_count = db_article_count - keyphrase_count
        self.target_and_nokeyphrase_count = target_count - \
                                            target_and_keyphrase_count
        self.notarget_and_keyphrase_count = keyphrase_count - \
                                            target_and_keyphrase_count
        self.notarget_and_nokeyphrase_count = self.nokeyphrase_count - \
                                            self.target_and_nokeyphrase_count


class KinderMinerResult(object):
    def __init__(self, target, p_value, target_and_keyphrase_ratio,
                                                    fet_p_value_and_ratio):
        self.target = target
        self.p_value = p_value
        self.target_and_keyphrase_ratio = target_and_keyphrase_ratio
        self.fet_p_value_and_ratio = fet_p_value_and_ratio


def compute_kinderminer_results(km_counts, p_cutoff=float('1e-5')):
    ret = list()
    for kmc in km_counts:
        # compute FET p
        targ_and_kp = kmc.target_and_keyphrase_count
        notarg_and_kp = kmc.notarget_and_keyphrase_count
        targ_and_nokp = kmc.target_and_nokeyphrase_count
        notarg_and_nokp = kmc.notarget_and_nokeyphrase_count
        odds_ratio,p_value = stats.fisher_exact([[targ_and_kp,notarg_and_kp], 
                                            [targ_and_nokp,notarg_and_nokp]], 
                                            alternative='greater')
        # eval
        denominator_sum = targ_and_kp + targ_and_nokp
        if p_value < p_cutoff and denominator_sum > 0:
            # compute ratio
            targ_and_kp_ratio = float(targ_and_kp) / float(denominator_sum)
            # computer FET p + ratio
            if p_value == 0.0:
                log_fet_p_value = -math.log10(float(1.0e-323))
            else:
                log_fet_p_value = -math.log10(float(p_value))
            log_ratio = math.log10(float(targ_and_kp_ratio))
            fet_p_value_and_ratio = log_fet_p_value + log_ratio
            ret.append(KinderMinerResult(kmc.target, 
                                         p_value, 
                                         targ_and_kp_ratio,
                                         fet_p_value_and_ratio))
    return ret

def main():
    COUNT_FILE = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]

    P_CUTOFF = float('1e-5')
    if len(sys.argv) > 3:
        P_CUTOFF = float(sys.argv[3])

    # read all the target counts with key phrase and alone
    counts_dict = dict()
    counts_list = list()
    with open(COUNT_FILE) as infile:
        header = infile.readline().strip().split('\t')
        targ_name_ind = header.index('target')
        targ_withkp_ind = header.index('target_with_keyphrase_count')
        targ_count_ind = header.index('target_count')
        kp_count_ind = header.index('keyphrase_count')
        total_count_ind = header.index('db_article_count')
        for line in infile:
            parts = line.strip().split('\t')
            name = parts[targ_name_ind]
            targ_tot = int(parts[targ_count_ind])
            kp_tot = int(parts[kp_count_ind])
            db_tot = int(parts[total_count_ind])
            targ_and_kp = int(parts[targ_withkp_ind])
            km_count = KinderMinerCount(name, 
                                        targ_and_kp, 
                                        targ_tot, 
                                        kp_tot, 
                                        db_tot)
            counts_dict[name] = km_count
            counts_list.append(km_count)

    # now compute the FET and ratios
    km_results = compute_kinderminer_results(counts_list, P_CUTOFF)
    
    # now sort by FET p-value and ratio
    km_results = sorted(km_results, 
                        key=lambda x: x.fet_p_value_and_ratio, 
                        reverse=True)
    
    # write to file
    with open(OUTPUT_FILE, 'w') as outfile:
        outfile.write('target\ttarget_with_keyphrase_count\t' +
                'target_count\tkeyphrase_count\tdb_article_count\t' +
                'fet_p_value\ttarget_and_keyphrase_ratio\t' + 
                'fet_p_value_and_ratio')
        outfile.write('\n')
        for kmr in km_results:
            targ = kmr.target
            targ_and_kp = counts_dict[targ].target_and_keyphrase_count
            targ_tot = counts_dict[targ].target_count
            kp_tot = counts_dict[targ].keyphrase_count
            db_tot = counts_dict[targ].db_article_count
            p_value = kmr.p_value
            tkp_ratio = kmr.target_and_keyphrase_ratio
            fet_ratio = kmr.fet_p_value_and_ratio
            outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(targ,
                                                            targ_and_kp,
                                                            targ_tot,
                                                            kp_tot,
                                                            db_tot,
                                                            p_value,
                                                            tkp_ratio,
                                                            fet_ratio))
            outfile.write('\n')


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()

    print(end - start)
