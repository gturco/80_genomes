from collections import Counter
from fuzzywuzzy import fuzz
import os
import re

#['hi my name is gene','hi my name is doug','hi my name is juli', 'hi my name is mikey']
####
### Cluster summary sheet
## CLUSTER #, Number of proteins in cluster, Number of genomes in the each CLUSTER, Avg number of times in
### geneome, % similarity between annoations, top five most common words


#### Cluster summary for each cluster
### Cluster Number
## Seen in 80 Genomes:
    ## Genome list
### Annoations are 70% identical
### Common Terms....


#Cluster_summary 
#    ##Number of geneomes cluster is in (have)
#    ##Average number of times seen in each geneome
#    ##Genenome distbution
#    % idenity of annoations (have)
#    listof most common words (have)
#    % of time seen, word

def write_cluster_summary(cluster_number,genome_cnt,anno_cnt,anno_similarity):
    """write summary of information for each cluster"""
    out_fh = open("cluster_sum/{0}_cluster_summary".format(cluster_number),"wb")
    total = sum(genome_cnt.values())
    avg_number = float(total)/(len(genome_cnt.keys()))
    genome_header = "seen in {0} of 80 genomes\n\n\n".format(len(genome_cnt.keys()))
    out_fh.write(genome_header)
    for genome in genome_cnt:
        genome_line = "{0}:{1}\n".format(genome,genome_cnt[genome])
        out_fh.write(genome_line)
    anno_header = "\n\n\nAnnotations are {0} % similar\n\n\n".format(anno_similarity)
    out_fh.write(anno_header)
    print anno_cnt
    top_five = set(anno_cnt)
    #top_five = ["{0}:{1}".format(key,value) for (key,value) in anno_cnt.most_common(3)]
    for anno in anno_cnt:
        anno_line = "{0}\n".format(anno)
        out_fh.write(anno_line)
    out_fh.close()
    return total, avg_number, ",".join(list(top_five))

def clean_anno(anno):
    new_anno = re.sub("'|\"","",anno)
    return new_anno

def count_anno(anno):
    """counts how many times the word shows in the cluster for each anno"""
    cnt = Counter()
    for sentence in anno:
        for word in sentence.split(' '):
            cnt[word] += 1
    ### sort dict by values
    return cnt

def percent_id(anno):
    """finds the average percent id for sentences in the cluster """
    all_sim = []
    anno.sort()
    new_annos = []
    for i,sentence in enumerate(anno):
        if i != len(anno)-1:
            simlarity = fuzz.token_sort_ratio(anno[i],anno[i+1])
            all_sim.append(simlarity)
            if simlarity <= 50:
                new_annos.append(sentence)
                new_annos.append(anno[i+1])
    if len(new_annos) == 0:
        new_annos.append(anno[1])
    avg = sum(all_sim)/len(all_sim)
    return avg, new_annos

def create_seq_dic(sequnce_file_dir):
    """creates the sequnce dic"""
    seq_dic = {}
    for sequnce_file in os.listdir(sequnce_file_dir):
        for line in open("{0}/{1}".format(sequnce_file_dir,sequnce_file)):
            if '>' in line:
                line_name = line.split(' ')[0]
                list_name = line_name[1:].strip().split('-')
                genome_id = '-'.join(list_name[:-1])
                i = int(list_name[-1])
                genome_key = "{0}_{1}".format(genome_id,str(i))
                seq_dic[genome_key] = ''
            elif len(line) > 1:
                seq_dic[genome_key] += line.strip()
            else: continue
    return seq_dic

def get_cluster_seq(seq_dic,genome_id,i,cluster_seq_out):
    """seq write outfile"""
    outfh = cluster_seq_out
    genome_key = "{0}_{1}".format(genome_id,i)
    seq = seq_dic[genome_key]
    line = ">{0}\n".format(genome_key)
    outfh.write(line)
    seq_line = "{0}\n".format(seq)
    outfh.write(seq_line)

def get_genome(gene_name):
    if 'scaffold' in gene_name:
        genome = ''.join(gene_name.split('_scaffold')[0])
    elif '_contig' in gene_name:
        genome = ''.join(gene_name.split('_contig')[0])
    else:
        genome = ''.join(gene_name.split('_gi|')[0])
    return genome

class cluster_line(object):
            
    def __init__(self, line):
        gene_name, i, anno = line.strip().split('\t')
        self.genome = get_genome(gene_name)
        self.clade_name = self.genome.split('_')[0]
        self.anno = clean_anno(anno)
        self.cluster_number = i


def main(cluster_dir,sequnce_file_dir):
    """creates a summary file for each cluster file """
    seq_dic = create_seq_dic(sequnce_file_dir)
    summary_file = open("cluster_sum_file_new_st.txt","wb")
    header = "##cluster_number\ttotal_in_cluster\tnumber_of_genomes\tavg_number\tsimilarity\ttop_three\n"
    summary_file.write(header)
    clade_file = open("clade_sum_st.txt","wb")
    all_clades = ['Halalkalicoccus','Haloarcula','Halobacterium','Halobiforma','Halococcus','Haloferax','Halogeometricum','Halomicrobium','Halopiger',
                'Haloquadratum','Halorhabdus','Halorubrum', 'Halosarcina' ,
                'Halosimplex', 'Haloterrigena', 'Halovivax', 'Natrialba', 'Natrinema',
                'Natronobacterium', 'Natronococcus', 'Natronolimnobius',
                'Natronomonas', 'Natronorubrum']
    header = 'cluster\tseen in genome\t{0}\n'.format('\t'.join(all_clades))
    clade_file.write(header)
    
    ### ALL NEW FiLES ####
    
    for cluster_file in os.listdir(cluster_dir):
        cluster_number = cluster_file.strip('_anno_cluster.txt')
        genome_cnt = Counter()
        clades = Counter()
        anno_list = []
        cluster_seq_out = open("protein_seq_{0}.fasta".format(cluster_number),'wb')
        for line in open("{0}/{1}".format(cluster_dir,cluster_file)):
            cluster = cluster_line(line)
            clades[cluster.clade_name] += 1
            genome_cnt[cluster.genome] += 1
            anno_list.append(cluster.anno)
            get_cluster_seq(seq_dic,cluster.genome,cluster.cluster_number,cluster_seq_out)
        cluster_seq_out.close()
        
        if len(anno_list) == 1: continue
        anno_similarity,anno_cnt = percent_id(anno_list)
        #anno_cnt = count_anno(anno_list)
        total, avg_number, top_five = write_cluster_summary(cluster_number,genome_cnt,anno_cnt,anno_similarity)
        summary_file_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(cluster_number,total,len(genome_cnt.keys()),avg_number,anno_similarity,top_five)
        summary_file.write(summary_file_line)
        
        clade_list = [cluster_number,str(len(genome_cnt.keys()))]
        for clade in all_clades:
            if clade in list(clades.keys()):
                col = clades[clade]
                clade_list.append(str(col))
            else:
                clade_list.append('0')
        print clade_list
        
        clade_file_line = "{0}\n".format('\t'.join(clade_list))
        clade_file.write(clade_file_line)



if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-q", dest="qfasta", help="path to genomic query fasta")




main('/Users/gt/Documents/code/eisen_fac/cluster_data/','/Users/gt/Dropbox/EightyHalophilesData/TranslationSequences_onlyCode/')






