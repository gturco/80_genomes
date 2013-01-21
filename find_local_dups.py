#Halosarcina_pallida_scaffold3_1|size567573_493863_494249
import os
from collections import defaultdict

def main(cluster_dir,pad):
    out_fh = open("cluster_dups.txt",'wb')
    for cluster_file in os.listdir(cluster_dir):
        cluster_number = cluster_file.split("_")[0]
        cluster_dic = parse_cluster_file("{0}/{1}".format(cluster_dir,cluster_file),pad)
        #print cluster_number
        cluster_total, cluster_overlap = count_clusters(cluster_dic)
        out_fh.write("{0}\t{1}\t{2}\n".format(cluster_number,cluster_overlap,cluster_total))
    out_fh.close()

def count_clusters(cluster_dic):
    cluster_total = 0
    cluster_overlap = 0
    for key in cluster_dic.keys():
        sub_cluster =  cluster_dic[key]
        if len(sub_cluster) == 1:
            cluster_total += 1
        else:
            cluster_total += len(sub_cluster)
            cluster_overlap += count_overlap(sub_cluster)
    return cluster_total, cluster_overlap



def find(x,y):
    """returns true if a set is not disjoint"""
    m = [x,y]
    m.sort()
    x,y = m
    return x[1] >= y[0] or x[1] >= y[1]

def count_overlap(genome_list):
    """counts number of overlaps for a list"""
    match_cnt = 0
    genome_list = [(int(s),int(e)) for (s,e) in genome_list]
    genome_list.sort()
    for gene in (genome_list):
        y = [gene] * len(genome_list)
        overlap = map(find,genome_list,y)
        if 1 < sum(overlap):
            match_cnt += 1
            #print overlap, genome_list,y
    return match_cnt

def parse_cluster_file(cluster_file,pad):
    """parse the cluster file"""
    cluster_dic = defaultdict(list)
    for line in open(cluster_file):
        line = line.strip().split('\t')[0]
        try:
            seqid =  "|".join([term for term in line.split('_') if '|' in term])
            gene_name, location = line.split(seqid)
            seqid = line.split("|")[0]
            start,end = map(int,location.split('_')[1:])
            key = "{0}_{1}".format(gene_name,seqid)
        except ValueError:
            gene_name = "_".join(line.split('_')[:-3])
            seqid,start,end = line.split('_')[-3:]
            start,end = map(int,[start,end])
            key = "{0}_{1}".format(gene_name,seqid)
        cluster_dic[key].append((max(0,start-pad),end+pad))
    return cluster_dic

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--cluster", dest="cluster_data", help="path to cluster_data_dir")
    parser.add_option("--pad", dest="padding", default=2000, help="bps to look up and downstream gene for local dups")
    (options, _) = parser.parse_args()


    main(options.cluster_data,options.padding)

#main("/Users/gt/Documents/code/eisen_fac/cluster_data/",2000)
#python find_local_dups --cluster "cluster_data/" --pad 2000

