from django.conf.urls.defaults import patterns, include, url
from Project.network.views import *
# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'cancerMiRNANetwork.views.home', name='home'),
    # url(r'^cancerMiRNANetwork/', include('cancerMiRNANetwork.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    ('^$',index),
    ('^index/$',index),
    ('^cancer/.*',cancer),
    ('.*/bp$',gene_ontology),
    ('.*/cc$',gene_ontology),
    ('.*/mf$',gene_ontology),
    ('.*/genes$',gene_listing),
    ('^overlapping_mirna$',overlapping_mirna),
    ('^overlapping_mirna_go$',overlapping_mirna_go),
    ('^common_mirna_and_go_term$',cytoscape_web_example),
    ('^cytoscape_web_common_mirna$',cytoscape_web_common_mirna),
    ('^cytoscape_web_common_mirna_go$',cytoscape_web_common_mirna_go),
    ('^inference/.*$',inference),
    ('^dataset/.*$',dataset),
    ('^cluster/.*$',cluster),
    ('^hallmark/.*$',hallmark),
    ('^hallmarks/.*$',hallmarks),
    ('^miRNA/.*$',miRNA),
    ('^count_em$',count_em),
    ('^count_em2$',count_em2),
    ('^specific_miRNA/.*$',specific_miRNA),
    ('^mirna_and_go_term$',mirna_and_go_term),
    ('^mirna_and_go_term_csv$',mirna_and_go_term_csv),
    ('^overlapping_mirna_go_csv$',overlapping_mirna_go_csv),
    ('^compendium/$',compendium),
    ('^firm/$',firm),
    ('^help/$',help_page),
    ('^citation/$',citation),
    ('^supTable8/$',supplementary_table_8_csv),
    ('^supTable9/$',supplementary_table_9_csv),
    ('^supTable10$',supplementary_table_10_csv),
    ('^overlap/.*$',overlap),
    ('^significant_overlapping_mirna_go_csv$',significant_overlapping_mirna_go_csv),
    ('^hallmarks_network_sif$',hallmarks_network_sif),
    ('^sig_overlap_network_sif$',sig_overlap_network_sif),
    ('^all_overlapping$',clusters_overlapping),
    ('^all_overlapping_hallmarks$',clusters_overlapping_hallmarks),
    ('^all_overlapping_hallmarks_mir2disease$',clusters_overlapping_hallmarks_mir2disease),
    ('^overlap_pita_targetscan$',clusters_overlapping_pita_targetscan),
    ('^overlap_pita_targetscan_cancer$',clusters_overlapping_pita_targetscan_cancer),
)

