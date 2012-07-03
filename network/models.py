from django.db import models

class Cancer(models.Model):
    name = models.CharField(max_length=200)
    short_name = models.CharField(max_length=100)
    publication = models.CharField(max_length=200)
    pmid = models.IntegerField()
    tissue = models.CharField(max_length=100)
    
    def __unicode__(self):
        return self.short_name

class Gene(models.Model):
    entrez_id = models.IntegerField()
    
    def __unicode__(self):
        return str(self.entrez_id)
    
class Coexpression_Cluster(models.Model):
    cancer = models.ForeignKey(Cancer)
    number = models.IntegerField()
    cluster_membership = models.ManyToManyField(Gene)
    
    def __unicode__(self):
        return self.cancer.__unicode__()+'.'+str(self.number)

class MiRNA(models.Model):
    name = models.CharField(max_length=100)
    mature_sequence_id = models.CharField(max_length=100)
    mature_sequence = models.CharField(max_length=100)
    pita_targets = models.ManyToManyField(Gene, related_name='MiRNAtoPITA')
    targetscan_targets = models.ManyToManyField(Gene,related_name='MiRNAtoTargetScan')
    
    def __unicode__(self):
        return self.name

class Validation(models.Model):
    cancer = models.ForeignKey(Cancer)
    mirna = models.ForeignKey(MiRNA)
    method = models.CharField(max_length=100)
    
    def __unicode__(self):
        return self.method+'.'+self.cancer.__unicode__()+'.'+self.mirna.__unicode__()

class Inferred_MiRNA(models.Model):
    coexpression_cluster = models.ForeignKey(Coexpression_Cluster)
    mirna = models.ForeignKey(MiRNA)
    method = models.CharField(max_length=100)
    target_genes = models.ManyToManyField(Gene, through='Target_Gene')
    
    def __unicode__(self):
        return self.coexpression_cluster.__unicode__()+'-'+self.mirna.__unicode__()

class Gene_Ontology(models.Model):
    go_id = models.CharField(max_length=100)
    term = models.CharField(max_length=500)
    category = models.CharField(max_length=100)
    annotated_genes = models.ManyToManyField(Gene)

    def __unicode__(self):
        return self.go_id

class Functional_Enrichment(models.Model):
    coexpression_cluster = models.ForeignKey(Coexpression_Cluster)
    gene_ontology = models.ForeignKey(Gene_Ontology)
    annotated_genes = models.ManyToManyField(Gene)

    def __unicode__(self):
        return self.coexpression_cluster.__unicode__()+'-'+self.gene_ontology.__unicode__()

class Target_Gene(models.Model):
    inferred_mirna = models.ForeignKey(Inferred_MiRNA)
    gene = models.ForeignKey(Gene)
    strand = models.CharField(max_length=10)
    site = models.CharField(max_length=10)
    start = models.IntegerField()
    match = models.FloatField()
    
    def __unicode__(self):
        return self.inferred_mirna.__unicode__()+' '+self.gene.__unicode__()

class Gene_Annotation(models.Model):
    gene = models.ForeignKey(Gene)
    category = models.CharField(max_length=100)
    annotation = models.CharField(max_length=1000)

    def __unicode__(self):
        return str(self.gene.entrez_id)+' = '+self.annotation

class Hallmarks_Of_Cancer(models.Model):
    coexpression_cluster = models.ForeignKey(Coexpression_Cluster)
    hallmark_of_cancer = models.CharField(max_length=1000)
    semantic_similarity = models.FloatField()
    
    def __unicode__(self):
        return self.coexpression_cluster.__unicode__()+' '+self.hallmark_of_cancer+'='+str(self.semantic_similarity)

class MiRNA_Target_GO_Term_Overlap(models.Model):
    coexpression_cluster = models.ForeignKey(Coexpression_Cluster)
    inferred_mirna = models.ForeignKey(Inferred_MiRNA)
    functional_enrichment = models.ForeignKey(Functional_Enrichment)
    p_value = models.FloatField()
    genes_overlapping = models.ManyToManyField(Gene)
    
    def __unicode__(self):
        return self.coexpression_cluster.__unicode__()+' '+self.inferred_mirna.__unicode__()+' '+self.gene_ontology.__unicode__()

