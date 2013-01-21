from collections import defaultdict
import re
import os
import sys

def clean_anno(anno):
    """cleans up anno so that you only get the product des this makes the data
    easier to work with"""
    m = re.findall("product=.*'|product=.*\n",anno)
    anno_new = m[0].strip().strip('product=')
    return anno_new

def create_cluster_dic(cluster_file):
    """create a dic of genome key id to cluster number"""
    fh = open(cluster_file)
    cluster = defaultdict(list)
    for line in fh:
        genome_id,scaffold,start,stop,cluster_number = line.strip().split('\t')
        accn = "{0}_{1}_{2}_{3}".format(genome_id, scaffold, start, stop)
        cluster[cluster_number].append(accn)
    return cluster

def create_anno_dic(genome_id, anno_file):
    """dic of gene_id to annoation for each genome"""
    anno_dic = {}
    i = 0
    for line in anno_file:
        scaffold,genebank,cds,start,stop,strand,rd,period,anno = line.split("\t")
        if cds != 'CDS': continue
        if 'genome_id' in anno: continue
        i += 1
        anno = clean_anno(anno)
        accn = "{0}_{1}_{2}_{3}".format(genome_id, scaffold, start, stop)
        anno_dic[accn] = (i,anno)
    return anno_dic

def join_dics(cluster,master_anno):
    """combines the anno dic with cluster dic keys are cluster numbers"""
    clust_anno = defaultdict(list)
    for cluster_number in cluster.keys():
        for gene in cluster[cluster_number]:
            gene_number,gene_anno = master_anno[gene]
            clust_anno[cluster_number].append((gene,gene_anno,gene_number))
    return clust_anno

class ClusterLine(object):
      
    def __init__(self, dirc, current_file):
        genome_fh = open("{0}/{1}".format(dirc,current_file))
        self.genome_id = current_file.split('.')[0]
        self.anno = create_anno_dic(self.genome_id, genome_fh)
        genome_fh.close()
 
def write_cluster_files(term,cluster,master_anno):
   """joins annno and clusters to write each cluster to a summary file"""
   clust_anno = join_dics(cluster,master_anno)
   for cluster_number in clust_anno.keys():
       anno_list = [a.upper() for (g,a,i) in clust_anno[cluster_number]]
       if term.upper() not in ','.join(anno_list): continue
       out_fh = open("cluster_data/{0}_anno_cluster.txt".format(cluster_number),'wb')
       for (gene, gene_anno, i) in clust_anno[cluster_number]:
           line = "{0}\t{1}\t{2}\n".format(gene,i,gene_anno)
           out_fh.write(line)
       out_fh.close()

   

def main(term,cluster_file,dirc):
    """opens annoation dir and cluster file and joins two files to dic """
    cluster = create_cluster_dic(cluster_file)
    master_anno = {}
    for current_file in os.listdir(dirc):
        print >>sys.stderr, "Opening .... {0}".format(current_file)
        c = ClusterLine(dirc,current_file)
        master_anno = dict(master_anno.items() + c.anno.items())
        ### adds annos to masteranno list for each genome
    #### write to file
    write_cluster_files(term,cluster,master_anno)


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("--term", dest="term", help="grab all annotations containing this term")
    parser.add_option("--cluster", dest="cluster_file", help="path to cluster file")
    parser.add_option("--anno", dest="anno_dir", help="directory containing annoation files")
    (options, _) = parser.parse_args()

    main(options.term, options.cluster_file, options.anno_dir)



#main('signal transduction','/Users/gt/Dropbox/Preliminary_ContextExplorer/EightyHalophiles/FilteredHomology/HomologyClusters_I20.txt','/Users/gt/Dropbox/Preliminary_ContextExplorer/EightyHalophiles/Annotations/')
#create_cluster_files.py --term 'signal transduction' --cluster
#'/Users/gt/Dropbox/Preliminary_ContextExplorer/EightyHalophiles/FilteredHomology/HomologyClusters_I20.txt',
#--anno '/Users/gt/Dropbox/Preliminary_ContextExplorer/EightyHalophiles/Annotations/')

