import portion as P
import itertools
from joblib import delayed,Parallel
import gffutils
import json
import argparse

def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def interval(t):
    return P.from_data([(True,i[0],i[1], True) for i in t])

def create_interval_dict_linear_time(gene,isoform_interval_dict):
    interval_set = set(isoform_interval_dict.keys())
    d = P.IntervalDict()
    union = P.empty()
    for transcript, inter in isoform_interval_dict.items():
        union = union | inter
    power_set_coords_dict = {}
    for p in P.iterate(union, step=1):
        s = list()
        for transcript, inter in isoform_interval_dict.items():
            if p in inter:
                s.append(transcript)
        s = repr(s)
        if s in power_set_coords_dict:
            power_set_coords_dict[s].append(p)
        else:
            power_set_coords_dict[s] = [p]
    for s, coords in power_set_coords_dict.items():
        d[interval(intervals_extract(coords))] = set(eval(s))
    return gene, d


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Write json file used for stitcher.py from a gtf file')
    parser.add_argument('-g','--gtf',metavar='gtf', type=str, help='Input gtf file')
    parser.add_argument('-d','--db', metavar='db', type=str, help='Intermediary database (db) file')
    parser.add_argument('-ji','--json_intervals', metavar='json', type=str, help='Output json file for coverage')
    parser.add_argument('-jr','--json_refskip', metavar='json', type=str, help='Output json file for refskip')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    gtffile = args.gtf
    dbfile = args.db
    jsonfile_1 = args.json_intervals
    jsonfile_2 = args.json_refskip
    threads = int(args.threads)
    print('Creating gtf database, this will take some time...')
    db = gffutils.create_db(gtffile, dbfile)
    isoform_interval_dict = {}
    isoform_refskip_dict = {}
    for gene in db.features_of_type('gene'):
        g_id = gene['gene_id'][0]
        isoform_interval_dict[g_id] = {}
        isoform_refskip_dict[g_id] = {}
        for transcript in db.children(gene, featuretype='transcript'):
            t_id = transcript['transcript_id'][0]
            isoform_interval_dict[g_id][t_id] = P.empty()
            isoform_refskip_dict[g_id][t_id] = P.empty()
            for exon in db.children(transcript, featuretype='exon'):
                isoform_interval_dict[g_id][t_id] = isoform_interval_dict[g_id][t_id] | P.closed(exon.start,exon.end)
            exon_list = [exon for exon in db.children(transcript, featuretype='exon')]
            if gene.strand == '+':
                for i in range(len(exon_list)-1):
                    isoform_refskip_dict[g_id][t_id] = isoform_refskip_dict[g_id][t_id] | P.closed(exon_list[i].end,exon_list[i+1].start)
            else:
                for i in range(len(exon_list)-1):
                    isoform_refskip_dict[g_id][t_id] = isoform_refskip_dict[g_id][t_id] | P.closed(exon_list[i+1].end,exon_list[i].start)
    print('Extracting unque isoform intervals')
    res = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(create_interval_dict_linear_time)(gene, transcript_intervals) for gene, transcript_intervals in isoform_interval_dict.items())
    isoform_unique_intervals = {k:v for k,v in res}
    isoform_unique_intervals_for_json_dump = {gene: {P.to_string(k):','.join(v) for k,v in d.items()} for gene,d in isoform_unique_intervals.items()}
    res_2 = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(create_interval_dict_linear_time)(gene, transcript_intervals) for gene, transcript_intervals in isoform_refskip_dict.items())
    isoform_unique_refskip = {k:v for k,v in res_2}
    isoform_unique_refskip_for_json_dump = {gene: {P.to_string(k):','.join(v) for k,v in d.items()} for gene,d in isoform_unique_refskip.items()}
    print('Writing unique isoform intervals to json file {}'.format(jsonfile_1))
    with open(jsonfile_1, 'w') as fp:
        json.dump(isoform_unique_intervals_for_json_dump, fp)
    print('Writing unique isoform refskip to json file {}'.format(jsonfile_2))
    with open(jsonfile_2, 'w') as fp:
        json.dump(isoform_unique_refskip_for_json_dump, fp)
