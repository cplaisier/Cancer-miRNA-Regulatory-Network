# Create your views here.
from django.http import HttpResponse
from django.template.loader import get_template
from django.template import Context
from django.db.models import Q
from network.models import *
from copy import deepcopy
import networkx as nx
from networkx.readwrite.graphml import write_graphml
from networkx.readwrite.gml import write_gml
import os, csv

# A method to write a sif and attributes file from a networkx graph
# node_attributes
def write_sif(G, node_attributes = [], edge_attributes = [], path='-'):
    # Write SIF file
    writeMe = []
    edges = G.edges()
    for edge in edges:
        writeMe.append(edge[0]+'\t'+G[edge[0]][edge[1]]['type']+'\t'+edge[1])
    sifFile = open(path+'.sif','w')
    sifFile.write('\n'.join(writeMe))
    sifFile.close()
    
    # Write node attribute files
    # Types = Double, Integer, String
    for attr in node_attributes:
        attribute = attr[0]
        attrType = attr[1]
        writeMe = []
        nodeNames = G.nodes()
        for node in nodeNames:
            writeMe.append(str(node)+' = '+str(G.node[node][attribute]))
        attrFile = open(path+'_'+str(attribute)+'_node.attr','w')
        attrFile.write(attribute+' (class='+attrType+')\n')
        attrFile.write('\n'.join(writeMe))
        attrFile.close()

    # Write edge attribute files
    # Types = Double, Integer, String
    for attr in edge_attributes:
        attribute = attr[0]
        attrType = attr[1]
        writeMe = []
        edges = G.edges()
        for edge in edges:
            if attribute in G[edge[0]][edge[1]]:
                writeMe.append(edge[0]+' ('+str(G[edge[0]][edge[1]][attribute])+') '+edge[1]+' = '+str(G[edge[0]][edge[1]][attribute]))
        attrFile = open(path+'_'+str(attribute)+'_edge.attr','w')
        attrFile.write(attribute+' (class='+attrType+')\n')
        attrFile.write('\n'.join(writeMe))
        attrFile.close()

def index(request):
    t = get_template('index.html')
    # Get all cancers and put them into an array of dictionaries
    hallmarks = {
            'self_sufficiency_in_growth_signals':'Self Sufficiency in Growth Signals',
            'insensitivity_to_antigrowth_signals':'Insensitivity to Antigrowth Signals',
            'evading_apoptosis':'Evading Apoptosis',
            'sustained_angiogenesis':'Sustained Angiogenesis',
            'tissue_invasion_and_metastasis':'Tissue Invasion and Meastasis',
            'reprogramming_energy_metabolism':'Reprogramming Energy Metabolism',
            'limitless_replicative_potential':'Limitless Replicative Potential',
            'genome_instability_and_mutation':'Genome Instability and Mutation',
            'evading_immune_detection':'Evading Immune Detection',
            'tumor_promoting_inflammation':'Tumor Promting Inflammation'
            }
    html = t.render(Context({'hallmarks':hallmarks}))
    return HttpResponse(html)

def cancer(request):
    t = get_template('cancer.html')
    cancer = request.path.split('/')[2]
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    im_cancer = Inferred_MiRNA.objects.filter(coexpression_cluster__cancer__tissue=cancer.lower())
    for im1 in im_cancer:
        clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
        if not clusterId in tmp:
            tmp[clusterId] = {}
            tmpUrls[clusterId] = {}
        if im1.method=='miRvestigator':
            if not 'miRvestigator' in tmp[clusterId]:
                tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
            else:
                tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
        if im1.method=='PITA':
            if not 'PITA' in tmp[clusterId]:
                tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
            else:
                tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
        if im1.method=='TargetScan':
            if not 'TargetScan' in tmp[clusterId]:
                tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
            else:
                tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

    fe_cancer = Functional_Enrichment.objects.filter(coexpression_cluster__cancer__tissue=cancer.lower())

    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    for fe1 in fe_cancer:
        clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
        if not clusterId in tmp:
            tmp[clusterId] = {}
        if fe1.gene_ontology.category=='biological_process':
            if not 'BP' in tmp[clusterId]:
                tmp[clusterId]['BP'] = []
            tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
        if fe1.gene_ontology.category=='molecular_function':
            if not 'MF' in tmp[clusterId]:
                tmp[clusterId]['MF'] = []
            tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
        if fe1.gene_ontology.category=='cellular_component':
            if not 'CC' in tmp[clusterId]:
                tmp[clusterId]['CC'] = []
            tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])

    # Create a wikipedia URL and pretty up the cancer name
    wikipedia = 'http://en.wikipedia.org/wiki/'
    if cancer=='Bladder' or cancer=='Brain' or cancer=='Breast' or cancer=='Colon' or cancer=='Head-Neck' or cancer=='Lung' or cancer=='Ovarian' or cancer=='Prostate' or cancer=='Renal':
        wikipedia = wikipedia + cancer + '_cancer'
        cancer = cancer + ' Cancer'
    elif cancer=='Pancreas':
        wikipedia = wikipedia + 'Pancreatic_cancer'
        cancer = 'Pancreatic Cancer'
    elif cancer=='Head-neck':
        wikipedia = wikipedia + 'Head_and_neck_cancer'
        cancer = 'Head and Neck Cancer'
    else:
        wikipedia = wikipedia + cancer

    html = t.render(Context({'cancer':cancer, 'wikipedia':wikipedia, 'entries':entries}))
    return HttpResponse(html)

def gene_ontology(request):
    t = get_template('list_go_terms.html')
    splitUp = request.path.split('/')
    cancer = splitUp[1].replace('_',' ')
    cluster = splitUp[2]
    termType = splitUp[3]
    if termType=='bp':
        termType = 'biological_process'
        category = 'Biological Process'
    elif termType=='cc':
        termType = 'cellular_component'
        category = 'Cellular Component'
    elif termType=='mf':
        termType = 'molecular_function'
        category = 'Molecular Function'
    fe_all = Functional_Enrichment.objects.filter(coexpression_cluster__number=int(cluster), coexpression_cluster__cancer__short_name=cancer, gene_ontology__category=termType)
    entries = []
    for fe1 in fe_all:
        num_genes = 0
        go1 = Gene_Ontology.objects.get(go_id=fe1.gene_ontology.go_id)
        go_genes = [i.entrez_id for i in go1.annotated_genes.all()]
        clust1 = Coexpression_Cluster.objects.get(cancer__short_name=cancer, number=int(cluster))
        clust_genes = [i.entrez_id for i in clust1.cluster_membership.all()]
        num_genes = len(set(go_genes).intersection(set(clust_genes)))
        if num_genes > 0:
            entries.append({ 'go_id':fe1.gene_ontology.go_id, 'term':fe1.gene_ontology.term, 'num_genes':num_genes, 'clust_genes':clust_genes, 'go_genes':go_genes })
    
    html = t.render(Context({ 'entries':entries, 'cancer':cancer, 'cluster':cluster, 'category': category }))
    return HttpResponse(html)

def gene_listing(request):
    t = get_template('go_term_list_genes.html')
    splitUp = request.path.split('/')
    cancer = splitUp[1]
    cluster = splitUp[2]
    term = splitUp[3]
    num_genes = 0
    go1 = Gene_Ontology.objects.get(go_id=term)
    go_genes = [i.entrez_id for i in go1.annotated_genes.all()]
    clust1 = Coexpression_Cluster.objects.get(cancer__short_name=cancer, number=int(cluster))
    clust_genes = [i.entrez_id for i in clust1.cluster_membership.all()]
    entries = []
    for gene in clust_genes:
        d1 = Gene_Annotation.objects.filter(gene__entrez_id=gene)
        cur = { 'entrez_id':gene }
        for i in d1:
            cur[i.category] = i.annotation
        if gene in go_genes:
            cur['yea_or_nea'] = 'Yes'
        else:
            cur['yea_or_nea'] = ''
        entries.append(cur)
        
    html = t.render(Context({ 'cancer': cancer, 'cluster': cluster, 'term': term, 'entries': entries }))
    return HttpResponse(html)

def overlapping_mirna(request):
    t = get_template('overlapping_mirna.html')
    # Grab all miRNA to co-expression cluster links
    im_all = Inferred_MiRNA.objects.all()
    mirna_counts = {}
    mirna_names = {}
    inferred_mirna = {}
    for im1 in im_all:
        if not im1.mirna.mature_sequence_id in mirna_counts:
            mirna_counts[im1.mirna.mature_sequence_id] = 1
            mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name.replace('hsa-','').replace('mir','miR')
            inferred_mirna[im1.mirna.mature_sequence_id] = [{'coexpression_cluster':im1.coexpression_cluster.__unicode__(), 'method':im1.method}]
        else:
            mirna_counts[im1.mirna.mature_sequence_id] += 1
            inferred_mirna[im1.mirna.mature_sequence_id].append({'coexpression_cluster':im1.coexpression_cluster.__unicode__(), 'method':im1.method})
    overlapping = {}
    for mirna in mirna_counts:
        if mirna_counts[mirna]>=2:
            overlapping[mirna_names[mirna]] = inferred_mirna[mirna]
    html = t.render(Context({ 'overlapping': overlapping }))
    return HttpResponse(html)

def overlapping_mirna_go(request):
    t = get_template('overlapping_mirna_go.html')
    # Grab all miRNA to co-expression cluster links
    im_all = Inferred_MiRNA.objects.all()
    mirna_counts = {}
    mirna_names = {}
    inferred_mirna = {}
    for im1 in im_all:
        if not im1.mirna.mature_sequence_id in mirna_counts:
            mirna_counts[im1.mirna.mature_sequence_id] = 0
            mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name
            inferred_mirna[im1.mirna.mature_sequence_id] = []
        mirna_counts[im1.mirna.mature_sequence_id] += 1
        inferred_mirna[im1.mirna.mature_sequence_id].append(deepcopy(im1))

    overlapping = []
    for mirna in mirna_counts:
        if mirna_counts[mirna]>=2:
            # Now check to see if there is any pairwise overlap between
            # any of the co-expression clusters with the same miRNA GO:<terms>
            functional_enrichments = {}
            for im1 in inferred_mirna[mirna]:
                fe_all = Functional_Enrichment.objects.filter(coexpression_cluster__number=im1.coexpression_cluster.number, coexpression_cluster__cancer__short_name=im1.coexpression_cluster.cancer, gene_ontology__category='biological_process')
                for fe1 in fe_all:
                    if not im1.coexpression_cluster.__unicode__() in functional_enrichments:
                        functional_enrichments[im1.coexpression_cluster.__unicode__()] = []
                    functional_enrichments[im1.coexpression_cluster.__unicode__()].append(fe1.gene_ontology.go_id)
            for fe1 in functional_enrichments.keys():
                functional_enrichments[fe1] = set(functional_enrichments[fe1])
            if len(functional_enrichments)>1:
                l1 = len(functional_enrichments)
                f1 = functional_enrichments.keys()
                for i1 in range(l1-1):
                    for i2 in range(i1+1,l1):
                        if not i1==i2:
                            o1 = functional_enrichments[f1[i1]].intersection(functional_enrichments[f1[i2]])
                            if len(o1):
                                terms = []
                                for go_id1 in o1:
                                    go1 = Gene_Ontology.objects.get(go_id=go_id1)
                                    terms.append({ 'go_id':go_id1, 'category':go1.category, 'term':go1.term })
                                overlapping.append({'mirna': mirna, 'mirna_name': mirna_names[mirna], 'gene_ontology':terms, 'cluster_1':f1[i1], 'cluster_2':f1[i2]})
    html = t.render(Context({ 'overlapping': overlapping }))
    return HttpResponse(html)

def cytoscape_web_example(request):
    t = get_template('cytoscape_web.html')
    html = t.render(Context({}))
    return HttpResponse(html)

def cytoscape_web_common_mirna(request):
    # Grab all miRNA to co-expression cluster links
    im_all = Inferred_MiRNA.objects.all()
    mirna_counts = {}
    mirna_names = {}
    inferred_mirna = {}
    for im1 in im_all:
        if not im1.mirna.mature_sequence_id in mirna_counts:
            mirna_counts[im1.mirna.mature_sequence_id] = 1
            mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name.replace('hsa-','').replace('mir','miR')
            inferred_mirna[im1.mirna.mature_sequence_id] = [{'cancer':im1.coexpression_cluster.cancer.tissue.capitalize(), 'cluster':im1.coexpression_cluster.__unicode__(), 'method':im1.method, 'mirna_name':im1.mirna.name}]
        else:
            mirna_counts[im1.mirna.mature_sequence_id] += 1
            inferred_mirna[im1.mirna.mature_sequence_id].append({'cancer':im1.coexpression_cluster.cancer.tissue.capitalize(), 'cluster':im1.coexpression_cluster.__unicode__(), 'method':im1.method, 'mirna_name':im1.mirna.name})
    overlapping = {}
    mirna_network = nx.Graph()
    for mirna in mirna_counts:
        if mirna_counts[mirna]>=10:
            for cluster in inferred_mirna[mirna]:
                #cluster = inferred_mirna[mirna][cluster]
                if not cluster['cancer'] in mirna_network.nodes():
                    mirna_network.add_node(cluster['cancer'],{'type':'tissue'})
                if not cluster['cluster'] in mirna_network.nodes():
                    mirna_network.add_node(cluster['cluster'], {'type':'cluster'})
                if not cluster['mirna_name'] in mirna_network.nodes():
                    mirna_network.add_node(cluster['mirna_name'], {'type':'mirna'})
                mirna_network.add_edge(cluster['cancer'], cluster['cluster'])
                mirna_network.add_edge(cluster['cluster'], cluster['mirna_name'])
    graphml = HttpResponse(content_type='application/xml')
    write_graphml(mirna_network, graphml)
    return graphml

def cytoscape_web_common_mirna_go(request):
    if not os.path.exists('/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer.sif'):
        # Grab all miRNA to co-expression cluster links
        im_all = Inferred_MiRNA.objects.all()
        mirna_counts = {}
        mirna_names = {}
        inferred_mirna = {}
        for im1 in im_all:
            if not im1.mirna.mature_sequence_id in mirna_counts:
                mirna_counts[im1.mirna.mature_sequence_id] = 0
                mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name
                inferred_mirna[im1.mirna.mature_sequence_id] = []
            mirna_counts[im1.mirna.mature_sequence_id] += 1
            inferred_mirna[im1.mirna.mature_sequence_id].append(deepcopy(im1))
        
        # Get all co-expression clusters with a hallmark of cancer (Semantic Similary >= 0.8)
        overlapping = []
        miRNAs = []
        coexpression_clusters = {}
        for mirna in mirna_counts:
            if mirna_counts[mirna]>=2:
                # Now check to see if there is any pairwise overlap between
                # any of the co-expression clusters with the same miRNA GO:<terms>
                functional_enrichments = {}
                for im1 in inferred_mirna[mirna]:
                    if not im1.coexpression_cluster.__unicode__() in coexpression_clusters:
                        coexpression_clusters[im1.coexpression_cluster.__unicode__()] = im1.coexpression_cluster.cancer.tissue.capitalize()
                    fe_all = Functional_Enrichment.objects.filter(coexpression_cluster__number=im1.coexpression_cluster.number, coexpression_cluster__cancer__short_name=im1.coexpression_cluster.cancer, gene_ontology__category='biological_process') | Functional_Enrichment.objects.filter(coexpression_cluster__number=im1.coexpression_cluster.number, coexpression_cluster__cancer__short_name=im1.coexpression_cluster.cancer, gene_ontology__category='molecular_function')
                    for fe1 in fe_all:
                        if not im1.coexpression_cluster.__unicode__() in functional_enrichments:
                            functional_enrichments[im1.coexpression_cluster.__unicode__()] = []
                        functional_enrichments[im1.coexpression_cluster.__unicode__()].append(fe1.gene_ontology.go_id)
                for fe1 in functional_enrichments.keys():
                    functional_enrichments[fe1] = set(functional_enrichments[fe1])
                if len(functional_enrichments)>1:
                    l1 = len(functional_enrichments)
                    f1 = functional_enrichments.keys()
                    for i1 in range(l1-1):
                        for i2 in range(i1+1,l1):
                            if not i1==i2:
                                o1 = functional_enrichments[f1[i1]].intersection(functional_enrichments[f1[i2]])
                                if len(o1):
                                    terms = []
                                    for go_id1 in o1:
                                        go1 = Gene_Ontology.objects.get(go_id=go_id1)
                                        terms.append({ 'go_id':go_id1, 'category':go1.category, 'term':go1.term })
                                    miRNAs.append(mirna)
                                    overlapping.append({'mirna': mirna, 'mirna_name': mirna_names[mirna], 'gene_ontology':terms, 'cluster_1':f1[i1], 'tissue_1':coexpression_clusters[f1[i1]], 'cluster_2':f1[i2], 'tissue_2':coexpression_clusters[f1[i2]]})

        network = nx.Graph()
        cc1 = {}
        for overlap in overlapping:
            if not overlap['tissue_1'] in network.nodes():
                network.add_node(overlap['tissue_1'], {'type':'tissue'})
            if not overlap['tissue_2'] in network.nodes():
                network.add_node(overlap['tissue_2'], {'type':'tissue'})
            if not overlap['cluster_1'] in network.nodes():
                network.add_node(overlap['cluster_1'], {'type':'cluster'})
            if not overlap['cluster_2'] in network.nodes():
                network.add_node(overlap['cluster_2'], {'type':'cluster'})
            network.add_edge(overlap['tissue_1'], overlap['cluster_1'], {'type':'tissue_cluster'} )
            network.add_edge(overlap['tissue_2'], overlap['cluster_2'], {'type':'tissue_cluster'})
            if not overlap['mirna_name'] in network.nodes():
                network.add_node(overlap['mirna_name'], {'type':'mirna', 'miRBase_id':overlap['mirna']})
            network.add_edge(overlap['cluster_1'], overlap['mirna_name'], {'type':'cluster_mirna'})
            network.add_edge(overlap['cluster_2'], overlap['mirna_name'], {'type':'cluster_mirna'})
            #for go1 in overlap['gene_ontology']:
            #    if not overlap['cluster_1'] in 
            #    go_ids = []
            #    terms = []
            #    categories = []
            #    if not go1['go_id'] in network.nodes():
            #        go_ids
            #        network.add_node(go1['go_id'], {'type':'go_id', 'category':go1['category'], 'term':go1['term']})
            #    network.add_edge(overlap['cluster_1'], go1['go_id'], {'type':'cluster_go_id'})
            #    network.add_edge(overlap['cluster_2'], go1['go_id'], {'type':'cluster_go_id'})
        
        # Hallmarks of cancer
        hoc_all = Hallmarks_Of_Cancer.objects.filter(semantic_similarity__gte=float(0.8))
        for hoc1 in hoc_all:
            im_all = Inferred_MiRNA.objects.filter(coexpression_cluster=hoc1.coexpression_cluster)
            if len(im_all) > 0:
                for im1 in im_all:
                    if im1.mirna.mature_sequence_id in miRNAs:
                        if not hoc1.hallmark_of_cancer in network.nodes():
                            network.add_node(hoc1.hallmark_of_cancer, {'type':'hallmark_of_cancer'})
                        if not [im1.mirna.name, hoc1.hallmark_of_cancer] in network.edges():
                            network.add_edge(im1.mirna.name, hoc1.hallmark_of_cancer, {'type':'mirna_hallmark_of_cancer'})
        
        # Write sif file
        write_sif(network, node_attributes=[['type','String']], edge_attributes=[['type','String']], path='/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer')

        # Build graphml from neworkx object
        outFile = open('/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer.gml','w')
        write_gml(network, outFile)
        outFile.close()
        gml = HttpResponse(content_type='application/xml')
        write_graphml(network, gml)
    else:
        gml = HttpResponse(content_type='application/xml')
        inFile = open('/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer.gml','r')
        [gml.write(line) for line in inFile.readlines() ] 
        inFile.close()
    return gml

def count_em(request):
    t = get_template('count_em.html')
    # Get all cancers and put them into an array of dictionaries
    fe_all = Functional_Enrichment.objects.filter(gene_ontology__category='biological_process')
    cc_all = []
    for fe1 in fe_all:
        cc1 = fe1.coexpression_cluster.__unicode__()
        if not cc1 in cc_all:
            cc_all.append(cc1)
    im_all = Inferred_MiRNA.objects.all()
    coincident_cc = []
    for im1 in im_all:
        cc1 = im1.coexpression_cluster.__unicode__()
        if cc1 in cc_all:
            if not cc1 in coincident_cc:
                coincident_cc.append(cc1)
    html = t.render(Context({'num_fe':len(cc_all), 'num_coincident':len(coincident_cc)}))
    return HttpResponse(html)

def count_em2(request):
    t = get_template('count_em2.html')
    all_validated = []
    # Get all miRvestigator inferred miRNAs that overlap with miR2Disease
    miRvestigator_validated = []
    im_all = Inferred_MiRNA.objects.filter(method='miRvestigator')
    for im1 in im_all:
        v1 = Validation.objects.filter(cancer__short_name=im1.coexpression_cluster.cancer.short_name, mirna__mature_sequence_id=im1.mirna.mature_sequence_id)
        if len(v1)>0 and not im1.coexpression_cluster.__unicode__() in miRvestigator_validated:
            miRvestigator_validated.append(im1.coexpression_cluster.__unicode__())
            if not im1.coexpression_cluster.__unicode__() in all_validated:
                all_validated.append(im1.coexpression_cluster.__unicode__())

    # Get all PITA inferred miRNAs that overlap with miR2Disease
    PITA_validated = []
    im_all = Inferred_MiRNA.objects.filter(method='PITA')
    for im1 in im_all:
        v1 = Validation.objects.filter(cancer__short_name=im1.coexpression_cluster.cancer.short_name, mirna__mature_sequence_id=im1.mirna.mature_sequence_id)
        if len(v1)>0 and not im1.coexpression_cluster.__unicode__() in PITA_validated:
            PITA_validated.append(im1.coexpression_cluster.__unicode__())
            if not im1.coexpression_cluster.__unicode__() in all_validated:
                all_validated.append(im1.coexpression_cluster.__unicode__())
    # Get all miRvestigator inferred miRNAs that overlap with miR2Disease
    TargetScan_validated = []
    im_all = Inferred_MiRNA.objects.filter(method='TargetScan')
    for im1 in im_all:
        v1 = Validation.objects.filter(cancer__short_name=im1.coexpression_cluster.cancer.short_name, mirna__mature_sequence_id=im1.mirna.mature_sequence_id)
        if len(v1)>0 and not im1.coexpression_cluster.__unicode__() in TargetScan_validated:
            TargetScan_validated.append(im1.coexpression_cluster.__unicode__())
            if not im1.coexpression_cluster.__unicode__() in all_validated:
                all_validated.append(im1.coexpression_cluster.__unicode__())
    html = t.render(Context({'miRvestigator':len(miRvestigator_validated), 'PITA':len(PITA_validated), 'TargetScan':len(TargetScan_validated), 'all':len(all_validated)}))
    return HttpResponse(html)

def inference(request):
    t = get_template('inference.html')
    splitUp = request.path.split('/')
    method = splitUp[2]
    # Get all cancers and put them into an array of dictionaries
    im_all = Inferred_MiRNA.objects.filter(method=method)
    miRNAs = []
    tmp = {}
    mirna2mat_seq_id = {}
    clusters = {}
    for im1 in im_all:
        cluster = im1.coexpression_cluster.__unicode__()
        mirna = im1.mirna.__unicode__()
        mat_seq_id = im1.mirna.mature_sequence_id
        if not cluster in tmp:
            tmp[cluster] = []
        if not mirna in tmp[cluster]:
            tmp[cluster].append(mirna)
            mirna2mat_seq_id[mirna] = mat_seq_id
        fe_dataset = Functional_Enrichment.objects.filter(coexpression_cluster=im1.coexpression_cluster)
        for fe1 in fe_dataset:
            if not cluster in clusters:
                clusters[cluster] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in clusters[cluster]:
                    clusters[cluster]['BP'] = []
                clusters[cluster]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in clusters[cluster]:
                    clusters[cluster]['MF'] = []
                clusters[cluster]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in clusters[cluster]:
                    clusters[cluster]['CC'] = []
                clusters[cluster]['CC'].append(str(fe1.gene_ontology.go_id))
    for cluster in tmp:
        splitUp = cluster.split('.')
        dataset = splitUp[0]
        number = splitUp[1]
        tmpDict = {'dataset':dataset, 'number':number, 'miRNA':dict(zip(tmp[cluster], [mirna2mat_seq_id[i] for i in tmp[cluster]])), 'url_bp':dataset.replace(' ','_')+'/'+str(number)+'/bp', 'url_mf':dataset.replace(' ','_')+'/'+str(number)+'/mf', 'url_cc':dataset.replace(' ','_')+'/'+str(number)+'/cc'}
        if cluster in clusters:
            if 'BP' in clusters[cluster]:
                tmpDict['BP'] = 'BP'
            if 'CC' in clusters[cluster]:
                tmpDict['CC'] = 'CC'
            if 'MF' in clusters[cluster]:
                tmpDict['MF'] = 'MF'
        miRNAs.append(deepcopy(tmpDict))
    html = t.render(Context({'method':method, 'entries':miRNAs}))
    return HttpResponse(html)

def dataset(request):
    t = get_template('dataset.html')
    splitUp = request.path.split('/')
    dataset = splitUp[2]
    # Get all cancers and put them into an array of dictionaries
    clusters = {}
    im_dataset = Inferred_MiRNA.objects.filter(coexpression_cluster__cancer__short_name=dataset)
    for im1 in im_dataset:
        cluster = im1.coexpression_cluster.number
        if not cluster in clusters:
            clusters[cluster] = {}
            clusters[cluster]['mirna'] = []
            clusters[cluster]['mat_seq_id'] = []
        clusters[cluster]['mirna'].append(im1.mirna.name)
        clusters[cluster]['mat_seq_id'].append(im1.mirna.mature_sequence_id)
    
    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    fe_dataset = Functional_Enrichment.objects.filter(coexpression_cluster__cancer__short_name=dataset)
    for fe1 in fe_dataset:
        cluster = fe1.coexpression_cluster.number
        if not cluster in clusters:
            clusters[cluster] = {}
        if fe1.gene_ontology.category=='biological_process':
            if not 'BP' in clusters[cluster]:
                clusters[cluster]['BP'] = []
            clusters[cluster]['BP'].append(str(fe1.gene_ontology.go_id))
        if fe1.gene_ontology.category=='molecular_function':
            if not 'MF' in clusters[cluster]:
                clusters[cluster]['MF'] = []
            clusters[cluster]['MF'].append(str(fe1.gene_ontology.go_id))
        if fe1.gene_ontology.category=='cellular_component':
            if not 'CC' in clusters[cluster]:
                clusters[cluster]['CC'] = []
            clusters[cluster]['CC'].append(str(fe1.gene_ontology.go_id))
    entries = []
    for entry in clusters:
        clusters[entry]['cluster'] = entry
        if 'mirna' in clusters[entry]:
            clusters[entry]['miRNAs'] = dict(zip(clusters[entry]['mirna'],clusters[entry]['mat_seq_id']))
        clusters[entry]['url_bp'] = dataset.replace(' ','_')+'/'+str(entry)+'/bp'
        clusters[entry]['url_mf'] = dataset.replace(' ','_')+'/'+str(entry)+'/mf'
        clusters[entry]['url_cc'] = dataset.replace(' ','_')+'/'+str(entry)+'/cc'
        entries.append(clusters[entry])
    c1 = Cancer.objects.get(short_name=dataset)
    cancer = c1.tissue
    wikipedia = 'http://en.wikipedia.org/wiki/'
    if cancer=='bladder' or cancer=='brain' or cancer=='breast' or cancer=='colon' or cancer=='lung' or cancer=='ovarian' or cancer=='prostate' or cancer=='renal':
        wikipedia = wikipedia + cancer + '_cancer'
        cancer = cancer + ' Cancer'
    elif cancer=='pancreas':
        wikipedia = wikipedia + 'Pancreatic_cancer'
        cancer = 'Pancreatic Cancer'
    elif cancer=='head-neck':
        wikipedia = wikipedia + 'Head_and_neck_cancer'
        cancer = 'Head and Neck Cancer'
    else:
        wikipedia = wikipedia
    c1.publication = c1.publication.strip('"')
    html = t.render(Context({'dataset':dataset, 'cancer':cancer, 'entries':entries, 'wikipedia':wikipedia, 'details':c1 }))
    return HttpResponse(html)

def cluster(request):
    t = get_template('cluster.html')
    splitUp = request.path.split('/')
    cancer = splitUp[2]
    cluster = splitUp[3]
    # Get all cancers and put them into an array of dictionaries
    c1 = Coexpression_Cluster.objects.get(cancer__short_name=cancer, number=int(cluster))
    clust_genes = [i.entrez_id for i in c1.cluster_membership.all()]
    # Regualtion
    miRNAs = {}
    regulation = {}
    im_all = Inferred_MiRNA.objects.filter(coexpression_cluster=c1)
    for im1 in im_all:
        if not im1.method in regulation:
            regulation[im1.method] = {}
        regulation[im1.method][im1.mirna.name] = im1.mirna.mature_sequence_id
    # Functional Enrichment
    go_terms = {}
    functional_enrichment = []
    fe_all = Functional_Enrichment.objects.filter(coexpression_cluster=c1)
    for fe1 in fe_all:
        num_genes = len(fe1.annotated_genes.all())
        if num_genes>0:
            functional_enrichment.append({ 'go_id':fe1.gene_ontology.go_id, 'go_term':fe1.gene_ontology.term, 'num_genes':len(fe1.annotated_genes.all()) })
    # Make overlap
    overlap = []
    o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=c1)
    for o1 in o_all:
        cluster_genes = c1.cluster_membership.all()
        if o1.inferred_mirna.method == 'miRvestigator':
            mirna_genes = set(o1.inferred_mirna.target_genes.all()).intersection(set(cluster_genes))
        elif o1.inferred_mirna.method == 'PITA':
            mirna_genes = set(o1.inferred_mirna.mirna.pita_targets.all()).intersection(set(cluster_genes))
        elif o1.inferred_mirna.method == 'TargetScan':
            mirna_genes = set(o1.inferred_mirna.mirna.targetscan_targets.all()).intersection(set(cluster_genes))
        overlap.append({'mat_seq_id':o1.inferred_mirna.mirna.mature_sequence_id, 'miRNA':o1.inferred_mirna.mirna.name, 'method':o1.inferred_mirna.method, 'miRNA_genes':len(mirna_genes), 'go_id':o1.functional_enrichment.gene_ontology.go_id, 'go_term':o1.functional_enrichment.gene_ontology.term, 'go_genes':len(o1.functional_enrichment.annotated_genes.all()), 'cluster_size':len(cluster_genes), 'overlap':len(o1.genes_overlapping.all()), 'p_value':o1.p_value })
    # Genes
    dataset = cancer
    c1 = Cancer.objects.get(short_name=cancer)
    cancer = c1.tissue
    wikipedia = 'http://en.wikipedia.org/wiki/'
    if cancer=='bladder' or cancer=='brain' or cancer=='breast' or cancer=='colon' or cancer=='lung' or cancer=='ovarian' or cancer=='prostate' or cancer=='renal':
        wikipedia = wikipedia + cancer + '_cancer'
        cancer = cancer + ' Cancer'
    elif cancer=='pancreas':
        wikipedia = wikipedia + 'Pancreatic_cancer'
        cancer = 'Pancreatic Cancer'
    elif cancer=='head-neck':
        wikipedia = wikipedia + 'Head_and_neck_cancer'
        cancer = 'Head and Neck Cancer'
    else:
        wikipedia = wikipedia
    c1.publication = c1.publication.strip('"')
    html = t.render(Context({'dataset':dataset, 'cancer':cancer, 'cluster':cluster, 'num_clust_genes': len(clust_genes), 'regulation':regulation, 'functional_enrichment':functional_enrichment, 'wikipedia':wikipedia, 'details':c1, 'overlaps':overlap })) #, 'genes':genes}))
    return HttpResponse(html)

def miRNA(request):
    t = get_template('miRNA_list_genes.html')
    splitUp = request.path.split('/')
    cancer = splitUp[2]
    cluster = int(splitUp[3])
    miRNA = splitUp[4]
    # Get miRvestigator target genes
    try:
 
        miRvestigator_target_genes = [i.entrez_id for i in Inferred_MiRNA.objects.get(coexpression_cluster__cancer__short_name=cancer, coexpression_cluster__number=cluster, mirna__mature_sequence_id=miRNA, method='miRvestigator').target_genes.all()]
    except:
        miRvestigator_target_genes = []
    # Get PITA target genes
    m1 = MiRNA.objects.get(mature_sequence_id=miRNA)
    PITA_target_genes = [i.entrez_id for i in m1.pita_targets.all()]
    # Get TargetScan target genes
    TargetScan_target_genes = [i.entrez_id for i in m1.targetscan_targets.all()]
    # Get cluster genes
    genes = [i.entrez_id for i in Coexpression_Cluster.objects.get(cancer__short_name=cancer, number=cluster).cluster_membership.all()]
    entries = []
    for gene in genes:
        d1 = Gene_Annotation.objects.filter(gene__entrez_id=gene)
        entry = { 'entrez_id':gene }
        for i in d1:
            entry[i.category] = i.annotation
        if gene in miRvestigator_target_genes:
            entry['target_miRvestigator'] = 'M'
        if gene in PITA_target_genes:
            entry['target_PITA'] = 'P'
        if gene in TargetScan_target_genes:
            entry['target_TargetScan'] = 'T'
        entries.append(entry)
    html = t.render(Context({'cancer':cancer, 'cluster':cluster, 'miRNA':m1.name, 'entries':entries}))
    return HttpResponse(html)

def specific_miRNA(request):
    t = get_template('specific_miRNA.html')
    splitUp = request.path.split('/')
    miRNA = splitUp[2]
    entries = []
    im_all = Inferred_MiRNA.objects.filter(mirna__mature_sequence_id=miRNA)
    for im1 in im_all:
        entries.append({'dataset':im1.coexpression_cluster.cancer.short_name,'cluster':im1.coexpression_cluster.number,'method':im1.method})
    html = t.render(Context({'miRNA':miRNA, 'entries':entries}))
    return HttpResponse(html)

def mirna_and_go_term(request):
    t = get_template('miRNA_and_GO_term.html')
    # Get all cancers and put them into an array of dictionaries
    fe_all = Functional_Enrichment.objects.filter(gene_ontology__category='biological_process')
    cc_all = []
    cc_all_dict = {}
    for fe1 in fe_all:
        cc1 = fe1.coexpression_cluster.__unicode__()
        if not cc1 in cc_all:
            cc_all.append(cc1)
            cc_all_dict[cc1] = {'funcEnrich':[]}
        cc_all_dict[cc1]['funcEnrich'].append(deepcopy(fe1))
    im_all = Inferred_MiRNA.objects.all()
    coincident_cc = []
    coincident_cc_dict = {}
    for im1 in im_all:
        cc1 = im1.coexpression_cluster.__unicode__()
        if cc1 in cc_all:
            if not cc1 in coincident_cc:
                coincident_cc.append(cc1)
                coincident_cc_dict[cc1] = {'funcEnrich':cc_all_dict[cc1]['funcEnrich'], 'miRNA':[]}
            coincident_cc_dict[cc1]['miRNA'].append(deepcopy(im1))
    entries = []
    for cc1 in coincident_cc_dict.keys():
        # dataset,cluster,miRNA,GO_id:GO_term
        entry = {'dataset':coincident_cc_dict[cc1]['funcEnrich'][0].coexpression_cluster.cancer.short_name, 'number':coincident_cc_dict[cc1]['funcEnrich'][0].coexpression_cluster.number, 'miRNAs':{}, 'funcEnrich':{}}
        # miRNAs -> miRNA name and mature_sequence_id
        for miRNA in coincident_cc_dict[cc1]['miRNA']:
            entry['miRNAs'][miRNA.mirna.name] = miRNA.mirna.mature_sequence_id
        # miRNAs -> miRNA name and mature_sequence_id
        for fe1 in coincident_cc_dict[cc1]['funcEnrich']:
            entry['funcEnrich'][fe1.gene_ontology.go_id] = { 'term':fe1.gene_ontology.term, 'category':fe1.gene_ontology.category }
        entries.append(entry)
    html = t.render(Context({ 'entries': entries }))
    return HttpResponse(html)

def mirna_and_go_term_csv(request):
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=coincident_go_terms.csv'
    # Get all cancers and put them into an array of dictionaries
    fe_all = Functional_Enrichment.objects.filter(gene_ontology__category='biological_process')
    cc_all = []
    cc_all_dict = {}
    for fe1 in fe_all:
        cc1 = fe1.coexpression_cluster.__unicode__()
        if not cc1 in cc_all:
            cc_all.append(cc1)
            cc_all_dict[cc1] = {'funcEnrich':[]}
        cc_all_dict[cc1]['funcEnrich'].append(deepcopy(fe1))
    im_all = Inferred_MiRNA.objects.all()
    coincident_cc = []
    coincident_cc_dict = {}
    for im1 in im_all:
        cc1 = im1.coexpression_cluster.__unicode__()
        if cc1 in cc_all:
            if not cc1 in coincident_cc:
                coincident_cc.append(cc1)
                coincident_cc_dict[cc1] = {'funcEnrich':cc_all_dict[cc1]['funcEnrich'], 'miRNA':[]}
            coincident_cc_dict[cc1]['miRNA'].append(deepcopy(im1))
    writer = csv.writer(response)
    entries = []
    writer.writerow(['Dataset','Cluster','GO:Terms','miRNAs','Validated'])
    for cc1 in coincident_cc_dict.keys():
        # dataset,cluster,miRNA,GO_id:GO_term
        entry = {'dataset':coincident_cc_dict[cc1]['funcEnrich'][0].coexpression_cluster.cancer.short_name, 'number':coincident_cc_dict[cc1]['funcEnrich'][0].coexpression_cluster.number, 'miRNAs':{}, 'funcEnrich':{}}
        # miRNAs -> miRNA name and mature_sequence_id
        validated = []
        for miRNA in coincident_cc_dict[cc1]['miRNA']:
            entry['miRNAs'][miRNA.mirna.name] = miRNA.mirna.mature_sequence_id
            v1 = Validation.objects.filter(cancer__short_name=entry['dataset'], mirna__mature_sequence_id=miRNA.mirna.mature_sequence_id)
            if not len(v1)==0:
                validated.append(miRNA.mirna.mature_sequence_id)
        # miRNAs -> miRNA name and mature_sequence_id
        for fe1 in coincident_cc_dict[cc1]['funcEnrich']:
            entry['funcEnrich'][fe1.gene_ontology.go_id] = { 'term':fe1.gene_ontology.term, 'category':fe1.gene_ontology.category }
        entries.append(entry)
        writer.writerow([entry['dataset'], entry['number'], '|'.join(entry['funcEnrich'].keys()), '|'.join(entry['miRNAs'].values()), '|'.join(validated)])
    return response

def compendium(request):
    t = get_template('compendium.html')
    html = t.render(Context())
    return HttpResponse(html)

def firm(request):
    t = get_template('firm.html')
    html = t.render(Context())
    return HttpResponse(html)

def help_page(request):
    t = get_template('help.html')
    html = t.render(Context())
    return HttpResponse(html)

def citation(request):
    t = get_template('citation.html')
    html = t.render(Context())
    return HttpResponse(html)

def hallmark(request):
    t = get_template('hallmark.html')
    hallmarkDict = {'evading_apoptosis':'Evading Apoptosis','evading_immune_detection':'Evading Immune Detection','genome_instability_and_mutation':'Genome Instability and Mutation','insensitivity_to_antigrowth_signals':'Insensitivity to Antigrowth Signals','limitless_replicative_potential':'Limitless Replicative Potential','reprogramming_energy_metabolism':'Reprogramming Energy Metabolism','self_sufficiency_in_growth_signals':'Self Sufficiency in Growth Signals','sustained_angiogenesis':'Sustained Angiogenesis','tissue_invasion_and_metastasis':'Tissue Invasion and Metastasis','tumor_promoting_inflammation':'Tumor Promting Inflammation'}
    hallmark = request.path.split('/')[2]
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    hmc_hallmark = Hallmarks_Of_Cancer.objects.filter(hallmark_of_cancer=hallmarkDict[hallmark], semantic_similarity__gte=0.8)
    for hmc1 in hmc_hallmark:
        im_hallmark = Inferred_MiRNA.objects.filter(coexpression_cluster=hmc1.coexpression_cluster)
        for im1 in im_hallmark:
            clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
                tmpUrls[clusterId] = {}
            if im1.method=='miRvestigator':
                if not 'miRvestigator' in tmp[clusterId]:
                    tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                else:
                    tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='PITA':
                if not 'PITA' in tmp[clusterId]:
                    tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='TargetScan':
                if not 'TargetScan' in tmp[clusterId]:
                    tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

        fe_hallmark = Functional_Enrichment.objects.filter(coexpression_cluster=hmc1.coexpression_cluster)

        # Get pathways for queried cancer
        # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
        for fe1 in fe_hallmark:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])

    # Create a wikipedia URL and pretty up the cancer name
    html = t.render(Context({'hallmark':hallmarkDict[hallmark], 'entries':entries}))
    return HttpResponse(html)

def hallmarks(request):
    t = get_template('hallmark.html')
    hallmarkDict = {'evading_apoptosis':'Evading Apoptosis','evading_immune_detection':'Evading Immune Detection','genome_instability_and_mutation':'Genome Instability and Mutation','insensitivity_to_antigrowth_signals':'Insensitivity to Antigrowth Signals','limitless_replicative_potential':'Limitless Replicative Potential','reprogramming_energy_metabolism':'Reprogramming Energy Metabolism','self_sufficiency_in_growth_signals':'Self Sufficiency in Growth Signals','sustained_angiogenesis':'Sustained Angiogenesis','tissue_invasion_and_metastasis':'Tissue Invasion and Metastasis','tumor_promoting_inflammation':'Tumor Promting Inflammation'}
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    hmc_hallmark = Hallmarks_Of_Cancer.objects.filter(semantic_similarity__gte=0.8)
    for hmc1 in hmc_hallmark:
        im_hallmark = Inferred_MiRNA.objects.filter(coexpression_cluster=hmc1.coexpression_cluster)
        for im1 in im_hallmark:
            clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
                tmpUrls[clusterId] = {}
            if im1.method=='miRvestigator':
                if not 'miRvestigator' in tmp[clusterId]:
                    tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                else:
                    tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='PITA':
                if not 'PITA' in tmp[clusterId]:
                    tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='TargetScan':
                if not 'TargetScan' in tmp[clusterId]:
                    tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

        fe_hallmark = Functional_Enrichment.objects.filter(coexpression_cluster=hmc1.coexpression_cluster)

        # Get pathways for queried cancer
        # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
        for fe1 in fe_hallmark:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])

    # Create a wikipedia URL and pretty up the cancer name
    html = t.render(Context({'hallmark':'all','entries':entries}))
    return HttpResponse(html)

def supplementary_table_8_csv(request):
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=coincident_miRNA_and_go_terms.csv'
    clusterDict = {}
    # Get all inferred miRNAs and check if they are validated by miR2Disease
    im_all = Inferred_MiRNA.objects.all()
    for im1 in im_all:
        if not im1.coexpression_cluster in clusterDict:
            clusterDict[im1.coexpression_cluster] = {}
        if not im1.method+'.miRNA' in clusterDict[im1.coexpression_cluster]:
            clusterDict[im1.coexpression_cluster][im1.method+'.miRNA'] = [im1.mirna.name]
            clusterDict[im1.coexpression_cluster][im1.method+'.mat_seq_id'] = [im1.mirna.mature_sequence_id]
        else:
            clusterDict[im1.coexpression_cluster][im1.method+'.miRNA'].append(im1.mirna.name)
            clusterDict[im1.coexpression_cluster][im1.method+'.mat_seq_id'].append(im1.mirna.mature_sequence_id)
        v_all = Validation.objects.filter(cancer__short_name=im1.coexpression_cluster.cancer.short_name, mirna__mature_sequence_id=im1.mirna.mature_sequence_id)
        if len(v_all) > 0:
            causal = 0
            for v1 in v_all:
                if v1.method=='miR2Disease:causal':
                    causal = 1
            if causal==1:
                clusterDict[im1.coexpression_cluster][im1.method+'.validated'] = 'causal'
            else:
                clusterDict[im1.coexpression_cluster][im1.method+'.validated'] = 'dysregulated'
    writer = csv.writer(response)
    writer.writerow(['Co-Expression Signature','miRvestigator miRNA','miRvestiagtor Validation','PITA miRNA','PITA Validation','TargetScan miRNA','TargetScan Validation'])
    for cluster in clusterDict:
        writeMe = [cluster]
        if 'miRvestigator.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['miRvestigator.miRNA']))
        else:
            writeMe.append('NA')
        if 'miRvestigator.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['miRvestigator.validated'])
        else:
            writeMe.append('NA')
        if 'PITA.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['PITA.miRNA']))
        else:
            writeMe.append('NA')
        if 'PITA.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['PITA.validated'])
        else:
            writeMe.append('NA')
        if 'TargetScan.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['TargetScan.miRNA']))
        else:
            writeMe.append('NA')
        if 'TargetScan.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['TargetScan.validated'])
        else:
            writeMe.append('NA')
        writer.writerow(writeMe)
    return response

## Co-incident miRNAs
def supplementary_table_9_csv(request):
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=coincident_miRNA_and_go_terms_sig_overlap_0.05.csv'
    clusterDict = {}
    # Get all inferred miRNAs and check if they are validated by miR2Disease
    im_all = Inferred_MiRNA.objects.all()
    for im1 in im_all:
        if not im1.coexpression_cluster in clusterDict:
            clusterDict[im1.coexpression_cluster] = {}
        if not im1.method+'.miRNA' in clusterDict[im1.coexpression_cluster]:
            clusterDict[im1.coexpression_cluster][im1.method+'.miRNA'] = [im1.mirna.name]
            clusterDict[im1.coexpression_cluster][im1.method+'.mat_seq_id'] = [im1.mirna.mature_sequence_id]
        else:
            clusterDict[im1.coexpression_cluster][im1.method+'.miRNA'].append(im1.mirna.name)
            clusterDict[im1.coexpression_cluster][im1.method+'.mat_seq_id'].append(im1.mirna.mature_sequence_id)
        v_all = Validation.objects.filter(cancer__short_name=im1.coexpression_cluster.cancer.short_name, mirna__mature_sequence_id=im1.mirna.mature_sequence_id)
        if len(v_all) > 0:
            causal = 0
            for v1 in v_all:
                if v1.method=='miR2Disease:causal':
                    causal = 1
            if causal==1:
                clusterDict[im1.coexpression_cluster][im1.method+'.validated'] = 'causal'
            else:
                clusterDict[im1.coexpression_cluster][im1.method+'.validated'] = 'dysregulated'
        mtgto_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=im1.coexpression_cluster, inferred_mirna=im1, p_value__lte=0.05)
        if len(mtgto_all)>0:
            clusterDict[im1.coexpression_cluster]['sig_overlap'] = []
            for mtgto1 in mtgto_all:
                clusterDict[im1.coexpression_cluster]['sig_overlap'].append(str(mtgto1.inferred_mirna.mirna.name)+'_'+str(mtgto1.functional_enrichment.gene_ontology.go_id))
    writer = csv.writer(response)
    writer.writerow(['Co-Expression Signature','miRvestigator miRNA','miRvestiagtor Validation','PITA miRNA','PITA Validation','TargetScan miRNA','TargetScan Validation','GO and miRNA Overlap'])
    for cluster in clusterDict:
        writeMe = [cluster]
        if 'miRvestigator.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['miRvestigator.miRNA']))
        else:
            writeMe.append('NA')
        if 'miRvestigator.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['miRvestigator.validated'])
        else:
            writeMe.append('NA')
        if 'PITA.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['PITA.miRNA']))
        else:
            writeMe.append('NA')
        if 'PITA.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['PITA.validated'])
        else:
            writeMe.append('NA')
        if 'TargetScan.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['TargetScan.miRNA']))
        else:
            writeMe.append('NA')
        if 'TargetScan.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['TargetScan.validated'])
        else:
            writeMe.append('NA')
        if 'sig_overlap' in clusterDict[cluster]:
            writeMe.append(';'.join(clusterDict[cluster]['sig_overlap']))
        else:
            writeMe.append('NA')
        writer.writerow(writeMe)
    return response

## Co-incident miRNAs
def supplementary_table_10_csv(request):
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=supplementary_table_10.csv'
    clusterDict = {}
    # Get all inferred miRNAs and check if they are validated by miR2Disease
    im_all = Inferred_MiRNA.objects.all()
    for im1 in im_all:
        if not im1.coexpression_cluster in clusterDict:
            clusterDict[im1.coexpression_cluster] = {}
        if not im1.method+'.miRNA' in clusterDict[im1.coexpression_cluster]:
            clusterDict[im1.coexpression_cluster][im1.method+'.miRNA'] = [im1.mirna.name]
            clusterDict[im1.coexpression_cluster][im1.method+'.mat_seq_id'] = [im1.mirna.mature_sequence_id]
        else:
            clusterDict[im1.coexpression_cluster][im1.method+'.miRNA'].append(im1.mirna.name)
            clusterDict[im1.coexpression_cluster][im1.method+'.mat_seq_id'].append(im1.mirna.mature_sequence_id)
        v_all = Validation.objects.filter(cancer__short_name=im1.coexpression_cluster.cancer.short_name, mirna__mature_sequence_id=im1.mirna.mature_sequence_id)
        if len(v_all) > 0:
            causal = 0
            for v1 in v_all:
                if v1.method=='miR2Disease:causal':
                    causal = 1
            if causal==1:
                clusterDict[im1.coexpression_cluster][im1.method+'.validated'] = 'causal'
            else:
                clusterDict[im1.coexpression_cluster][im1.method+'.validated'] = 'dysregulated'
        mtgto_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=im1.coexpression_cluster, inferred_mirna=im1, p_value__lte=0.05)
        if len(mtgto_all)>0:
            clusterDict[im1.coexpression_cluster]['sig_overlap'] = []
            for mtgto1 in mtgto_all:
                clusterDict[im1.coexpression_cluster]['sig_overlap'].append(str(mtgto1.inferred_mirna.mirna.name)+'_'+str(mtgto1.functional_enrichment.gene_ontology.go_id))
        hmc_all = Hallmarks_Of_Cancer.objects.filter(coexpression_cluster=im1.coexpression_cluster, hallmark_of_cancer='Tissue Invasion and Metastasis', semantic_similarity__gte=0.8)
        if len(hmc_all)>0:
            clusterDict[im1.coexpression_cluster]['metastatic'] = 1
    writer = csv.writer(response)
    writer.writerow(['Co-Expression Signature','miRvestigator miRNA','miRvestiagtor Validation','PITA miRNA','PITA Validation','TargetScan miRNA','TargetScan Validation','GO and miRNA Overlap','Metastatic'])
    for cluster in clusterDict:
        writeMe = [cluster]
        if 'miRvestigator.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['miRvestigator.miRNA']))
        else:
            writeMe.append('NA')
        if 'miRvestigator.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['miRvestigator.validated'])
        else:
            writeMe.append('NA')
        if 'PITA.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['PITA.miRNA']))
        else:
            writeMe.append('NA')
        if 'PITA.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['PITA.validated'])
        else:
            writeMe.append('NA')
        if 'TargetScan.miRNA' in clusterDict[cluster]:
            writeMe.append('_'.join(clusterDict[cluster]['TargetScan.miRNA']))
        else:
            writeMe.append('NA')
        if 'TargetScan.validated' in clusterDict[cluster]:
            writeMe.append(clusterDict[cluster]['TargetScan.validated'])
        else:
            writeMe.append('NA')
        if 'sig_overlap' in clusterDict[cluster]:
            writeMe.append(';'.join(clusterDict[cluster]['sig_overlap']))
        else:
            writeMe.append('NA')
        if 'metastatic' in clusterDict[cluster]:
            writeMe.append('Y')
        else:
            writeMe.append('N')
        writer.writerow(writeMe)
    return response

def overlap(request):
    # http://cmrn.systemsbiology.net/overlap/AD%20Lung%20Bhattacharjee/41/MIMAT0000690/TargetScan/GO:0016477
    t = get_template('overlap.html')
    splitUp = request.path.split('/')
    cancer = splitUp[2]
    cluster = int(splitUp[3])
    miRNA = splitUp[4]
    method = splitUp[5]
    go_id = splitUp[6]
    go_term = Gene_Ontology.objects.get(go_id=go_id).term
    # Get miRvestigator target genes
    try:
         miRvestigator_target_genes = [i.entrez_id for i in Inferred_MiRNA.objects.get(coexpression_cluster__cancer__short_name=cancer, coexpression_cluster__number=cluster, mirna__mature_sequence_id=miRNA, method='miRvestigator').target_genes.all()]
    except:
        miRvestigator_target_genes = []
    # Get PITA target genes
    m1 = MiRNA.objects.get(mature_sequence_id=miRNA)
    PITA_target_genes = [i.entrez_id for i in m1.pita_targets.all()]
    # Get TargetScan target genes
    TargetScan_target_genes = [i.entrez_id for i in m1.targetscan_targets.all()]
    # Get cluster genes
    c1 = Coexpression_Cluster.objects.get(cancer__short_name=cancer, number=cluster)
    genes = [i.entrez_id for i in c1.cluster_membership.all()]
    go_genes = [i.entrez_id for i in Functional_Enrichment.objects.get(coexpression_cluster=c1, gene_ontology__go_id=go_id).annotated_genes.all()]
    entries = []
    for gene in genes:
        d1 = Gene_Annotation.objects.filter(gene__entrez_id=gene)
        entry = { 'entrez_id':gene }
        for i in d1:
            entry[i.category] = i.annotation
        if gene in miRvestigator_target_genes:
            entry['target_miRvestigator'] = 'M'
        if gene in PITA_target_genes:
            entry['target_PITA'] = 'P'
        if gene in TargetScan_target_genes:
            entry['target_TargetScan'] = 'T'
        if gene in go_genes:
            entry['go_gene'] = 'Yes'
        else:
            entry['go_gene'] = ''
        entries.append(entry)
    html = t.render(Context({'cancer':cancer, 'cluster':cluster, 'miRNA':m1.name, 'entries':entries, 'method':method, 'go_id':go_id, 'go_term':go_term}))
    return HttpResponse(html)

def significant_overlapping_mirna_go_csv(request):
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=sig_0.05_overlaping_coincident_miRNA_and_go_terms.csv'
    # Grab all miRNA to co-expression cluster links
    im_all = Inferred_MiRNA.objects.all()
    mirna_counts = {}
    mirna_names = {}
    inferred_mirna = {}
    for im1 in im_all:
        if not im1.mirna.mature_sequence_id in mirna_counts:
            mirna_counts[im1.mirna.mature_sequence_id] = 0
            mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name
            inferred_mirna[im1.mirna.mature_sequence_id] = []
        mirna_counts[im1.mirna.mature_sequence_id] += 1
        inferred_mirna[im1.mirna.mature_sequence_id].append(deepcopy(im1))

    overlapping = []
    for mirna in mirna_counts:
        if mirna_counts[mirna]>=2:
            # Now check to see if there is any pairwise overlap between
            # any of the co-expression clusters with the same miRNA GO:<terms>
            functional_enrichments = {}
            for im1 in inferred_mirna[mirna]:
                o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=im1.coexpression_cluster, inferred_mirna=im1, p_value__lte=0.05)
                for o1 in o_all:
                    if len(o1.genes_overlapping.all()) > 0:
                        fe1 = o1.functional_enrichment
                        if not im1.coexpression_cluster.__unicode__() in functional_enrichments:
                            functional_enrichments[im1.coexpression_cluster.__unicode__()] = []
                        functional_enrichments[im1.coexpression_cluster.__unicode__()].append(fe1.gene_ontology.go_id)
            for fe1 in functional_enrichments.keys():
                functional_enrichments[fe1] = set(functional_enrichments[fe1])
            if len(functional_enrichments)>1:
                l1 = len(functional_enrichments)
                f1 = functional_enrichments.keys()
                for i1 in range(l1-1):
                    for i2 in range(i1+1,l1):
                        if not i1==i2:
                            o1 = functional_enrichments[f1[i1]].intersection(functional_enrichments[f1[i2]])
                            if len(o1):
                                terms = []
                                for go_id1 in o1:
                                    go1 = Gene_Ontology.objects.get(go_id=go_id1)
                                    terms.append({ 'go_id':go_id1, 'category':go1.category, 'term':go1.term })
                                overlapping.append({'mirna': mirna, 'mirna_name': mirna_names[mirna], 'gene_ontology':terms, 'cluster_1':f1[i1], 'cluster_2':f1[i2]})
    
    writer = csv.writer(response)
    writer.writerow(['Dataset.1','Cluster.1','Dataset.2','Cluster.2','GO:Terms','miRNA','miRBase ID'])
    for overlap in overlapping:
        c1 = overlap['cluster_1'].split('.')
        c2 = overlap['cluster_2'].split('.')
        writer.writerow([c1[0],c1[1],c2[0],c2[1],'|'.join([i['go_id'] for i in overlap['gene_ontology']]),overlap['mirna'],overlap['mirna_name']])
    return response

def overlapping_mirna_go_csv(request):
    response = HttpResponse(mimetype='text/csv')
    response['Content-Disposition'] = 'attachment; filename=sig_overlaping_coincident_miRNA_and_go_terms.csv'
    # Grab all miRNA to co-expression cluster links
    im_all = Inferred_MiRNA.objects.all()
    mirna_counts = {}
    mirna_names = {}
    inferred_mirna = {}
    for im1 in im_all:
        if not im1.mirna.mature_sequence_id in mirna_counts:
            mirna_counts[im1.mirna.mature_sequence_id] = 0
            mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name
            inferred_mirna[im1.mirna.mature_sequence_id] = []
        mirna_counts[im1.mirna.mature_sequence_id] += 1
        inferred_mirna[im1.mirna.mature_sequence_id].append(deepcopy(im1))

    overlapping = []
    for mirna in mirna_counts:
        if mirna_counts[mirna]>=2:
            # Now check to see if there is any pairwise overlap between
            # any of the co-expression clusters with the same miRNA GO:<terms>
            functional_enrichments = {}
            for im1 in inferred_mirna[mirna]:
                fe_all = Functional_Enrichment.objects.filter(coexpression_cluster=im1.coexpression_cluster, gene_ontology__category='biological_process')
                for fe1 in fe_all:
                    if len(fe1.annotated_genes.all())>0:
                        if not im1.coexpression_cluster.__unicode__() in functional_enrichments:
                            functional_enrichments[im1.coexpression_cluster.__unicode__()] = []
                        functional_enrichments[im1.coexpression_cluster.__unicode__()].append(fe1.gene_ontology.go_id)
            for fe1 in functional_enrichments.keys():
                functional_enrichments[fe1] = set(functional_enrichments[fe1])
            if len(functional_enrichments)>1:
                l1 = len(functional_enrichments)
                f1 = functional_enrichments.keys()
                for i1 in range(l1-1):
                    for i2 in range(i1+1,l1):
                        if not i1==i2:
                            o1 = functional_enrichments[f1[i1]].intersection(functional_enrichments[f1[i2]])
                            if len(o1):
                                terms = []
                                for go_id1 in o1:
                                    go1 = Gene_Ontology.objects.get(go_id=go_id1)
                                    terms.append({ 'go_id':go_id1, 'category':go1.category, 'term':go1.term })
                                overlapping.append({'mirna': mirna, 'mirna_name': mirna_names[mirna], 'gene_ontology':terms, 'cluster_1':f1[i1], 'cluster_2':f1[i2]})
    
    writer = csv.writer(response)
    writer.writerow(['Dataset.1','Cluster.1','Dataset.2','Cluster.2','GO:Terms','miRNA','miRBase ID'])
    for overlap in overlapping:
        c1 = overlap['cluster_1'].split('.')
        c2 = overlap['cluster_2'].split('.')
        writer.writerow([c1[0],c1[1],c2[0],c2[1],'|'.join([i['go_id'] for i in overlap['gene_ontology']]),overlap['mirna'],overlap['mirna_name']])
    return response

def hallmarks_network_sif(request):
    if not os.path.exists('/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer_V2.sif'):
        # Grab all miRNA to co-expression cluster links
        im_all = Inferred_MiRNA.objects.all()
        mirna_counts = {}
        mirna_names = {}
        inferred_mirna = {}
        for im1 in im_all:
            if not im1.mirna.mature_sequence_id in mirna_counts:
                mirna_counts[im1.mirna.mature_sequence_id] = 0
                mirna_names[im1.mirna.mature_sequence_id] = im1.mirna.name
                inferred_mirna[im1.mirna.mature_sequence_id] = []
            mirna_counts[im1.mirna.mature_sequence_id] += 1
            inferred_mirna[im1.mirna.mature_sequence_id].append(deepcopy(im1))
        
        # Get all co-expression clusters with a hallmark of cancer (Semantic Similary >= 0.8)
        overlapping = []
        miRNAs = []
        coexpression_clusters = {}
        for mirna in mirna_counts:
            if mirna_counts[mirna]>=2:
                # Now check to see if there is any pairwise overlap between
                # any of the co-expression clusters with the same miRNA GO:<terms>
                functional_enrichments = {}
                for im1 in inferred_mirna[mirna]:
                    if not im1.coexpression_cluster.__unicode__() in coexpression_clusters:
                        coexpression_clusters[im1.coexpression_cluster.__unicode__()] = im1.coexpression_cluster.cancer.tissue.capitalize()
                    
                    o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=im1.coexpression_cluster, inferred_mirna=im1, p_value__lte=0.05)
                    for o1 in o_all:
                        if len(o1.genes_overlapping.all()) > 0:
                            fe1 = o1.functional_enrichment
                            if not im1.coexpression_cluster.__unicode__() in functional_enrichments:
                                functional_enrichments[im1.coexpression_cluster.__unicode__()] = []
                            functional_enrichments[im1.coexpression_cluster.__unicode__()].append(fe1.gene_ontology.go_id)
                for fe1 in functional_enrichments.keys():
                    functional_enrichments[fe1] = set(functional_enrichments[fe1])
                if len(functional_enrichments)>1:
                    l1 = len(functional_enrichments)
                    f1 = functional_enrichments.keys()
                    for i1 in range(l1-1):
                        for i2 in range(i1+1,l1):
                            if not i1==i2:
                                o1 = functional_enrichments[f1[i1]].intersection(functional_enrichments[f1[i2]])
                                if len(o1):
                                    terms = []
                                    for go_id1 in o1:
                                        go1 = Gene_Ontology.objects.get(go_id=go_id1)
                                        terms.append({ 'go_id':go_id1, 'category':go1.category, 'term':go1.term })
                                    miRNAs.append(mirna)
                                    overlapping.append({'mirna': mirna, 'mirna_name': mirna_names[mirna], 'gene_ontology':terms, 'cluster_1':f1[i1], 'tissue_1':coexpression_clusters[f1[i1]], 'cluster_2':f1[i2], 'tissue_2':coexpression_clusters[f1[i2]]})

        network = nx.Graph()
        cc1 = {}
        for overlap in overlapping:
            if not overlap['tissue_1'] in network.nodes():
                network.add_node(overlap['tissue_1'], {'type':'tissue'})
            if not overlap['tissue_2'] in network.nodes():
                network.add_node(overlap['tissue_2'], {'type':'tissue'})
            if not overlap['cluster_1'] in network.nodes():
                network.add_node(overlap['cluster_1'], {'type':'cluster'})
            if not overlap['cluster_2'] in network.nodes():
                network.add_node(overlap['cluster_2'], {'type':'cluster'})
            network.add_edge(overlap['tissue_1'], overlap['cluster_1'], {'type':'tissue_cluster'} )
            network.add_edge(overlap['tissue_2'], overlap['cluster_2'], {'type':'tissue_cluster'})
            if not overlap['mirna_name'] in network.nodes():
                network.add_node(overlap['mirna_name'], {'type':'mirna', 'miRBase_id':overlap['mirna']})
            network.add_edge(overlap['cluster_1'], overlap['mirna_name'], {'type':'cluster_mirna'})
            network.add_edge(overlap['cluster_2'], overlap['mirna_name'], {'type':'cluster_mirna'})
            #for go1 in overlap['gene_ontology']:
            #    if not overlap['cluster_1'] in 
            #    go_ids = []
            #    terms = []
            #    categories = []
            #    if not go1['go_id'] in network.nodes():
            #        go_ids
            #        network.add_node(go1['go_id'], {'type':'go_id', 'category':go1['category'], 'term':go1['term']})
            #    network.add_edge(overlap['cluster_1'], go1['go_id'], {'type':'cluster_go_id'})
            #    network.add_edge(overlap['cluster_2'], go1['go_id'], {'type':'cluster_go_id'})
        
        # Hallmarks of cancer
        hoc_all = Hallmarks_Of_Cancer.objects.filter(semantic_similarity__gte=float(0.8))
        for hoc1 in hoc_all:
            im_all = Inferred_MiRNA.objects.filter(coexpression_cluster=hoc1.coexpression_cluster)
            if len(im_all) > 0:
                for im1 in im_all:
                    if im1.mirna.mature_sequence_id in miRNAs:
                        if not hoc1.hallmark_of_cancer in network.nodes():
                            network.add_node(hoc1.hallmark_of_cancer, {'type':'hallmark_of_cancer'})
                        if not [im1.mirna.name, hoc1.hallmark_of_cancer] in network.edges():
                            network.add_edge(im1.mirna.name, hoc1.hallmark_of_cancer, {'type':'mirna_hallmark_of_cancer'})
        
        # Write sif file
        write_sif(network, node_attributes=[['type','String']], edge_attributes=[['type','String']], path='/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer_V2')

        # Build graphml from neworkx object
        outFile = open('/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer_V2.gml','w')
        write_gml(network, outFile)
        outFile.close()
        gml = HttpResponse(content_type='application/xml')
        write_graphml(network, gml)
    else:
        gml = HttpResponse(content_type='application/xml')
        inFile = open('/opt/bitnami/apps/django/django_projects/Project/network/static/hallmark_of_cancer_V2.gml','r')
        [gml.write(line) for line in inFile.readlines() ] 
        inFile.close()
    return gml

def clusters_overlapping(request):
    t = get_template('overlapping.html')
    # Get all clusters where overlap is significant between miRNA and GO term
    o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(p_value__lte=0.05)
    overlap = []
    clusters = {}
    for o1 in o_all:
        tmpName = o1.coexpression_cluster.cancer.short_name+' '+str(o1.coexpression_cluster.number)
        g_all = o1.genes_overlapping.all()
        if not tmpName in clusters and not len(g_all)==0:
            clusters[tmpName] = o1.coexpression_cluster
        tmpName2 = o1.coexpression_cluster.cancer.short_name+' '+str(o1.coexpression_cluster.number)+' '+o1.inferred_mirna.mirna.mature_sequence_id+' '+o1.inferred_mirna.method
        if not len(g_all)==0 and not tmpName2 in overlap:
            overlap.append(tmpName2)
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    for cluster in clusters:
        im_cancer = Inferred_MiRNA.objects.filter(coexpression_cluster=clusters[cluster])
        for im1 in im_cancer:
            # Only if signficant overlap exists
            o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=clusters[cluster], inferred_mirna=im1, p_value__lte=0.05)
            if len(o_all)>0:
                tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
                if tmpName in clusters:
                    clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
                    tmpName2 = clusters[cluster].cancer.short_name+' '+str(clusters[cluster].number)+' '+im1.mirna.mature_sequence_id+' '+im1.method
                    if not clusterId in tmp:
                        tmp[clusterId] = {}
                        tmpUrls[clusterId] = {}
                    if im1.method=='miRvestigator' and tmpName2 in overlap:
                        if not 'miRvestigator' in tmp[clusterId]:
                            tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                        else:
                            tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
                    if im1.method=='PITA' and tmpName2 in overlap:
                        if not 'PITA' in tmp[clusterId]:
                            tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                        else:
                            tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
                    if im1.method=='TargetScan' and tmpName2 in overlap:
                        if not 'TargetScan' in tmp[clusterId]:
                            tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                        else:
                            tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

    fe_cancer = Functional_Enrichment.objects.all()

    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    for fe1 in fe_cancer:
        tmpName = fe1.coexpression_cluster.cancer.short_name+' '+str(fe1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    miRNAs = 0
    miRNAs_matSeqIds = []
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            for matSeqID in tmp[entry]['miRvestigator']:
                if not matSeqID in miRNAs_matSeqIds:
                    miRNAs += 1
                    miRNAs_matSeqIds.append(matSeqID)
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            for matSeqID in tmp[entry]['PITA']:
                if not matSeqID in miRNAs_matSeqIds:
                    miRNAs += 1
                    miRNAs_matSeqIds.append(matSeqID)
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            for matSeqID in tmp[entry]['TargetScan']:
                if not matSeqID in miRNAs_matSeqIds:
                    miRNAs += 1
                    miRNAs_matSeqIds.append(matSeqID)
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])
    html = t.render(Context({'entries':entries,'miRNAs':miRNAs}))
    return HttpResponse(html)

def clusters_overlapping_hallmarks(request):
    t = get_template('overlapping.html')
    # Get all clusters where overlap is significant between miRNA and GO term
    o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(p_value__lte=0.05)
    clusters = {}
    for o1 in o_all:
        tmpName = o1.coexpression_cluster.cancer.short_name+' '+str(o1.coexpression_cluster.number)
        g_all = o1.genes_overlapping.all()
        hc_all = Hallmarks_Of_Cancer.objects.filter(coexpression_cluster=o1.coexpression_cluster, semantic_similarity__gte=0.8)
        if not tmpName in clusters and not len(g_all)==0 and not len(hc_all)==0:
            clusters[tmpName] = o1.coexpression_cluster
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    for cluster in clusters:
        im_cancer = Inferred_MiRNA.objects.filter(coexpression_cluster=clusters[cluster])
        for im1 in im_cancer:
            # Only if signficant overlap exists
            o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=clusters[cluster], inferred_mirna=im1, p_value__lte=0.05)
            if len(o_all)>0:
                tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
                if tmpName in clusters:
                    clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
                    if not clusterId in tmp:
                        tmp[clusterId] = {}
                        tmpUrls[clusterId] = {}
                    if im1.method=='miRvestigator':
                        if not 'miRvestigator' in tmp[clusterId]:
                            tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                        else:
                            tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
                    if im1.method=='PITA':
                        if not 'PITA' in tmp[clusterId]:
                            tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                        else:
                            tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
                    if im1.method=='TargetScan':
                        if not 'TargetScan' in tmp[clusterId]:
                            tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                        else:
                            tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

    fe_cancer = Functional_Enrichment.objects.all()

    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    for fe1 in fe_cancer:
        tmpName = fe1.coexpression_cluster.cancer.short_name+' '+str(fe1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    miRNAs = 0
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            miRNAs += 1
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            miRNAs += 1
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            miRNAs += 1
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])
    html = t.render(Context({'entries':entries,'miRNAs':miRNAs}))
    return HttpResponse(html)

def clusters_overlapping_hallmarks_mir2disease(request):
    t = get_template('overlapping.html')
    # Get all clusters where overlap is significant between miRNA and GO term
    o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(p_value__lte=0.05)
    clusters = {}
    for o1 in o_all:
        tmpName = o1.coexpression_cluster.cancer.short_name+' '+str(o1.coexpression_cluster.number)
        g_all = o1.genes_overlapping.all()
        hc_all = Hallmarks_Of_Cancer.objects.filter(coexpression_cluster=o1.coexpression_cluster, semantic_similarity__gte=0.8)
        v_all = Validation.objects.filter(cancer=o1.coexpression_cluster.cancer,mirna=o1.inferred_mirna.mirna)
        if not tmpName in clusters and not len(g_all)==0 and not len(hc_all)==0 and not len(v_all)==0:
            clusters[tmpName] = o1.coexpression_cluster
           
    # Get Inferred miRNAs 
    tmp = {}
    tmpUrls = {}
    for cluster in clusters:
        im_cancer = Inferred_MiRNA.objects.filter(coexpression_cluster=clusters[cluster])
        for im1 in im_cancer:
            # Only if signficant overlap exists
            o_all = MiRNA_Target_GO_Term_Overlap.objects.filter(coexpression_cluster=clusters[cluster], inferred_mirna=im1, p_value__lte=0.05)
            if len(o_all)>0:
                tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
                if tmpName in clusters:
                    clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
                    if not clusterId in tmp:
                        tmp[clusterId] = {}
                        tmpUrls[clusterId] = {}
                    if im1.method=='miRvestigator':
                        if not 'miRvestigator' in tmp[clusterId]:
                            tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                        else:
                            tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
                    if im1.method=='PITA':
                        if not 'PITA' in tmp[clusterId]:
                            tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                        else:
                            tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
                    if im1.method=='TargetScan':
                        if not 'TargetScan' in tmp[clusterId]:
                            tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                        else:
                            tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

    fe_cancer = Functional_Enrichment.objects.all()

    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    for fe1 in fe_cancer:
        tmpName = fe1.coexpression_cluster.cancer.short_name+' '+str(fe1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    miRNAs = 0
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            miRNAs += 1
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            miRNAs += 1
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            miRNAs += 1
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])
    html = t.render(Context({'entries':entries,'miRNAs':miRNAs}))
    return HttpResponse(html)

def clusters_overlapping_pita_targetscan(request):
    t = get_template('overlapping.html')
    im_all_pita = Inferred_MiRNA.objects.filter(method='PITA')
    clusters = {}
    for im1 in im_all_pita:
        if len(Inferred_MiRNA.objects.filter(method='TargetScan', coexpression_cluster=im1.coexpression_cluster, mirna=im1.mirna))>=1:
            tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
            clusters[tmpName] = im1.coexpression_cluster 
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    im_cancer = Inferred_MiRNA.objects.all()
    for im1 in im_cancer:
        tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
                tmpUrls[clusterId] = {}
            if im1.method=='miRvestigator':
                if not 'miRvestigator' in tmp[clusterId]:
                    tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                else:
                    tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='PITA':
                if not 'PITA' in tmp[clusterId]:
                    tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='TargetScan':
                if not 'TargetScan' in tmp[clusterId]:
                    tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

    fe_cancer = Functional_Enrichment.objects.all()

    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    for fe1 in fe_cancer:
        tmpName = fe1.coexpression_cluster.cancer.short_name+' '+str(fe1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])
    html = t.render(Context({'entries':entries}))
    return HttpResponse(html)

def clusters_overlapping_pita_targetscan_cancer(request):
    t = get_template('overlapping.html')
    im_all_pita = Inferred_MiRNA.objects.filter(method='PITA')
    clusters = {}
    for im1 in im_all_pita:
        im_all = Inferred_MiRNA.objects.filter(method='TargetScan', coexpression_cluster__cancer=im1.coexpression_cluster.cancer, mirna=im1.mirna)
        if len(im_all)>=1:
            tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
            clusters[tmpName] = im1.coexpression_cluster
            for im2 in im_all:
                tmpName = im2.coexpression_cluster.cancer.short_name+' '+str(im2.coexpression_cluster.number)
                if not tmpName in clusters:
                    clusters[tmpName] = im2.coexpression_cluster
            
    # Get Inferred miRNA where cancer queried cancer
    tmp = {}
    tmpUrls = {}
    im_cancer = Inferred_MiRNA.objects.all()
    for im1 in im_cancer:
        tmpName = im1.coexpression_cluster.cancer.short_name+' '+str(im1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = im1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
                tmpUrls[clusterId] = {}
            if im1.method=='miRvestigator':
                if not 'miRvestigator' in tmp[clusterId]:
                    tmp[clusterId]['miRvestigator'] = {im1.mirna.mature_sequence_id: { 'name': im1.mirna.name.lstrip('hsa-').replace('mir','miR') }}
                else:
                    tmp[clusterId]['miRvestigator'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='PITA':
                if not 'PITA' in tmp[clusterId]:
                    tmp[clusterId]['PITA'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['PITA'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}
            if im1.method=='TargetScan':
                if not 'TargetScan' in tmp[clusterId]:
                    tmp[clusterId]['TargetScan'] = {im1.mirna.mature_sequence_id:{'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}}
                else:
                    tmp[clusterId]['TargetScan'][im1.mirna.mature_sequence_id] = {'name':im1.mirna.name.lstrip('hsa-').replace('mir','miR')}

    fe_cancer = Functional_Enrichment.objects.all()

    # Get pathways for queried cancer
    # A bad use-case exists where the link will be available even if there are no genes in the cluster that map to the go_id!
    for fe1 in fe_cancer:
        tmpName = fe1.coexpression_cluster.cancer.short_name+' '+str(fe1.coexpression_cluster.number)
        if tmpName in clusters:
            clusterId = fe1.coexpression_cluster.__unicode__().replace('_',' ').replace('.',' - ')
            if not clusterId in tmp:
                tmp[clusterId] = {}
            if fe1.gene_ontology.category=='biological_process':
                if not 'BP' in tmp[clusterId]:
                    tmp[clusterId]['BP'] = []
                tmp[clusterId]['BP'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='molecular_function':
                if not 'MF' in tmp[clusterId]:
                    tmp[clusterId]['MF'] = []
                tmp[clusterId]['MF'].append(str(fe1.gene_ontology.go_id))
            if fe1.gene_ontology.category=='cellular_component':
                if not 'CC' in tmp[clusterId]:
                    tmp[clusterId]['CC'] = []
                tmp[clusterId]['CC'].append(str(fe1.gene_ontology.go_id))

    entries = []
    for entry in tmp:
        tmp[entry]['cluster'] = entry
        splitUp = entry.split(' - ')
        tmp[entry]['dataset'] = splitUp[0]
        tmp[entry]['number'] = splitUp[1]
        if entry in tmp and 'miRvestigator' in tmp[entry].keys():
            for matSeqID in tmp[entry]['miRvestigator']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['miRvestigator'][matSeqID]['validation'] = 1
        if entry in tmp and 'PITA' in tmp[entry].keys():
            for matSeqID in tmp[entry]['PITA']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['PITA'][matSeqID]['validation'] = 1
        if entry in tmp and 'TargetScan' in tmp[entry].keys():
            for matSeqID in tmp[entry]['TargetScan']:
                v1 = Validation.objects.filter(cancer__short_name=tmp[entry]['dataset'], mirna__mature_sequence_id=matSeqID)
                if not len(v1)==0:
                    tmp[entry]['TargetScan'][matSeqID]['validation'] = 1
        tmp[entry]['url_bp'] = entry.replace(' - ','/').replace(' ','_')+'/bp'
        tmp[entry]['url_mf'] = entry.replace(' - ','/').replace(' ','_')+'/mf'
        tmp[entry]['url_cc'] = entry.replace(' - ','/').replace(' ','_')+'/cc'
        entries.append(tmp[entry])
    html = t.render(Context({'entries':entries}))
    return HttpResponse(html)

def sig_overlap_network_sif(request):
    if not os.path.exists('/opt/bitnami/apps/django/django_projects/Project/network/static/sig_overlap_network_0.05.gml'):
        ## Output a network of coexpression clusters that have significant number of miRNA target genes that are also GO annotated
        coexpression_clusters = {}
        mtgto_all = MiRNA_Target_GO_Term_Overlap.objects.filter(p_value__lte=0.05)
        for mtgto1 in mtgto_all:
            # Exclude if significant with zero genes overlapping
            og_all = mtgto1.genes_overlapping.all()
            # Include if co-expression cluster enriched for hallmark of cancer
            im1 = mtgto1.inferred_mirna
            hoc_all = Hallmarks_Of_Cancer.objects.filter(semantic_similarity__gte=float(0.8), coexpression_cluster=im1.coexpression_cluster)
            if len(og_all)>0 and len(hoc_all)>0:
                if not mtgto1.coexpression_cluster in coexpression_clusters:
                    coexpression_clusters[mtgto1.coexpression_cluster] = []
                if not im1.mirna in coexpression_clusters[mtgto1.coexpression_cluster]:
                    coexpression_clusters[mtgto1.coexpression_cluster].append(im1.mirna)

        # Now build the network
        network = nx.Graph()
        for cluster in coexpression_clusters:
            # Add tissue node if not already present
            tissue = cluster.cancer.tissue.capitalize()
            if not tissue in network.nodes():
                network.add_node(tissue, {'type':'tissue'})
            # Add cluster node if not already present
            clusterName = cluster.__unicode__()
            if not clusterName in network.nodes():
                network.add_node(clusterName, {'type':'cluster'})
            # Add edge between tissue and cluster
            if not (tissue, clusterName) in network.edges():
                network.add_edge(tissue, clusterName, {'type':'tissue_cluster'} )
            mirnaName = '_'.join([m1.name for m1 in coexpression_clusters[cluster]])
            mirbase = '_'.join([m1.mature_sequence_id for m1 in coexpression_clusters[cluster]])
            if not mirnaName in network.nodes():
                network.add_node(mirnaName, {'type':'mirna', 'miRBase_id':mirbase})
            if not (clusterName, mirnaName) in network.edges():
                validated = 0
                for mirna in coexpression_clusters[cluster]:
                    validated += len(Validation.objects.filter(cancer=cluster.cancer, mirna=mirna))
                if validated > 0:
                    network.add_edge(clusterName, mirnaName, {'type':'cluster_mirna','validated':0})
                else:
                    network.add_edge(clusterName, mirnaName, {'type':'cluster_mirna','validated':1})
            # Hallmarks of cancer
            hoc_all = Hallmarks_Of_Cancer.objects.filter(semantic_similarity__gte=float(0.8), coexpression_cluster=cluster)
            for hoc1 in hoc_all:
                if not hoc1.hallmark_of_cancer in network.nodes():
                    network.add_node(hoc1.hallmark_of_cancer, {'type':'hallmark_of_cancer'})
                if not (mirnaName, hoc1.hallmark_of_cancer) in network.edges():
                    network.add_edge(mirnaName, hoc1.hallmark_of_cancer, {'type':'mirna_hallmark_of_cancer'})
        
        # Write sif file
        write_sif(network, node_attributes=[['type','String']], edge_attributes=[['type','String'],['validated','String']], path='/opt/bitnami/apps/django/django_projects/Project/network/static/sig_overlap_network_0.05')

        # Build graphml from neworkx object
        outFile = open('/opt/bitnami/apps/django/django_projects/Project/network/static/sig_overlap_network_0.05.gml','w')
        write_gml(network, outFile)
        outFile.close()
        gml = HttpResponse(content_type='application/xml')
        write_graphml(network, gml)
    else:
        gml = HttpResponse(content_type='application/xml')
        inFile = open('/opt/bitnami/apps/django/django_projects/Project/network/static/sig_overlap_network_0.05.gml','r')
        [gml.write(line) for line in inFile.readlines() ] 
        inFile.close()
    return gml

