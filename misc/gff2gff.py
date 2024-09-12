#!/usr/bin/env python3
#
# This is script was contributed by Philip Ashton (https://github.com/flashton2003)
# with the intention to convert from the many GFF flavours to the Ensembl GFF3 format
# supported by bcftools csq. Please expand as necessary and send a pull request.
# 
# See also
#   https://github.com/samtools/bcftools/issues/530#issuecomment-268278248
#   https://github.com/samtools/bcftools/issues/1078#issuecomment-527484831
#   https://github.com/samtools/bcftools/issues/1208#issuecomment-620642373
#
import sys
import gffutils

class PseudoFeature():
    def __init__(self, chrom, start, stop, strand):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand
        

class FeatureGroup():
    def __init__(self, gene_id):
        ## see here https://github.com/samtools/bcftools/issues/530#issuecomment-614410603
        self.id = gene_id
        ## the below should be lists of instances of Feature class.
        self.gene = None # one feature group should have a single gene
        self.transcript = None
        ## assuming here that ncRNA is more like transcript than exon
        ## i.e. there is only one ncRNA feature per feature group.
        self.ncRNA = None
        self.exons = []
        self.CDSs = []

    def add_gene(self, feature):
        self.gene = feature

    def add_transcript(self, feature):
        self.transcript = feature # one feature group should have a single gene

    def add_exon(self, feature):
        self.exons.append(feature)

    def add_cds(self, feature):
        self.CDSs.append(feature)

    def add_ncRNA(self, feature):
        self.ncRNA = feature


def get_args():
    if len(sys.argv) != 3:
        print('Usage: gff2gff.py <gff_inhandle> </path/where/you/want/to/save/gffutils_db>')
        sys.exit()
    else:
        return sys.argv[1], sys.argv[2]

def group_features(db):
    feature_groups = {}
    for feature in db.all_features():
        if feature.featuretype not in ('repeat_region', 'regulatory', 'stem_loop', 'gene_component_region', 'repeat_region'):
            if feature.featuretype == 'gene':
                ## for any particular feature group, the gene doesn't necassarily come
                ## before the CDS in the gff, so need to have the if/else logic
                gene_id = feature.id.split('.')[0]
                if gene_id not in feature_groups:
                    feature_groups[gene_id] = FeatureGroup(gene_id)
                    feature_groups[gene_id].add_gene(feature)
                else:
                    feature_groups[gene_id].add_gene(feature)
            ## in this gff, the transcripts are referred to as mRNA
            ## not all genes have mRNA records
            elif feature.featuretype == 'mRNA':
                gene_id = feature.id.split('.')[0]
                # print(feature.id)
                # print(gene_id)
                feature_groups[gene_id].add_transcript(feature)
            elif feature.featuretype == 'exon':
                gene_id = feature.attributes['locus_tag'][0]
                if gene_id not in feature_groups:
                    feature_groups[gene_id] = FeatureGroup(gene_id)
                    feature_groups[gene_id].add_exon(feature)
                else:
                    feature_groups[gene_id].add_exon(feature)

            elif feature.featuretype == 'CDS':
                gene_id = feature.attributes['locus_tag'][0]
                if gene_id not in feature_groups:
                    feature_groups[gene_id] = FeatureGroup(gene_id)
                    feature_groups[gene_id].add_cds(feature)
                else:
                    feature_groups[gene_id].add_cds(feature)

            elif feature.featuretype == 'ncRNA':
                gene_id = feature.attributes['locus_tag'][0]
                # print(gene_id)
                if gene_id not in feature_groups:
                    feature_groups[gene_id] = FeatureGroup(gene_id)
                    feature_groups[gene_id].add_ncRNA(feature)
                else:
                    feature_groups[gene_id].add_ncRNA(feature)
    return feature_groups

def make_transcript_pseudofeature(CDSs):
    ## takes list of CDSs, list len will often be 1, but will sometimes be 2
    start = min([x.start for x in CDSs])
    stop = max([x.stop for x in CDSs])
    chrom = CDSs[0].chrom
    ## check that both on the same strand
    assert len(set([x.strand for x in CDSs])) == 1
    strand = CDSs[0].strand
    pf = PseudoFeature(chrom, start, stop, strand)
    return pf


def check_feature_groups(feature_groups):
    ## checking that feature group has a gene, at least one CDS, and a transcript
    ## if doesn't have a transcript (as many don't), then make a pseudo feature corresponding to a transcript
    ## only checking that have transcript for everything
    for gene_id in feature_groups:
        ## don't want to analyse ncRNAs here, but need to include them up 
        ## to this point so that we know that that gene is an ncRNA gene
        if feature_groups[gene_id].ncRNA != None:
            continue
        assert feature_groups[gene_id].gene != None
        assert len(feature_groups[gene_id].CDSs) > 0
        if feature_groups[gene_id].transcript == None:
            ## using pseudofeature because the gffutils feature says it's not really intended
            ## for direct instantiation by users
            feature_groups[gene_id].transcript = make_transcript_pseudofeature(feature_groups[gene_id].CDSs)
       
def print_features(feature_groups):
    for gene_id in feature_groups:
        feature_group = feature_groups[gene_id]
        if feature_group.ncRNA != None:
            continue
        print('###')
        name = feature_group.gene.attributes['Name'][0]
        gene_attributes = f'ID=gene:{gene_id};Name={name};biotype=protein_coding;gene_id:{gene_id}'
        print(feature_group.gene.chrom, 'EMBL', 'gene', feature_group.gene.start, feature_group.gene.stop, '.', feature_group.gene.strand, '.', gene_attributes, sep = '\t')
        
        transcript_attributes = f'ID=transcript:{gene_id};Parent=gene:{gene_id};Name={name};biotype=protein_coding;transcript_id={gene_id}'
        print(feature_group.transcript.chrom, 'EMBL', 'transcript', feature_group.transcript.start, feature_group.transcript.stop, '.', feature_group.transcript.strand, '.', transcript_attributes, sep = '\t')
        
        cds_attributes =  f'Parent=transcript:{gene_id};Name={name}'

        for c in feature_group.CDSs:
            print(c.chrom, 'EMBL', 'CDS', c.start, c.stop, '.', c.strand, '0', cds_attributes, sep = '\t')

def main():
    '''
    Ensembl gff should have one gene and one transcript per "feature group"
    Then, can have multiple CDS/exons
    read in from the input gff, one FeatureGroup instance has one gene, one transcript and (potentially)
    multiple CDS/exon
    Exons aren't printed, as not needed bt bcftools csq
    '''
    gff_handle, db_handle = get_args()
    fn = gffutils.example_filename(gff_handle)
    db = gffutils.create_db(fn, dbfn=db_handle, force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
    feature_groups = group_features(db)
    check_feature_groups(feature_groups)
    print_features(feature_groups)



## gff input is generated from a genbank file by bp_genbank2gff3.pl

if __name__ == '__main__':
    main()

