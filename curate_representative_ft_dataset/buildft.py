import click
import logging
from rich.logging import RichHandler
from Bio import Phylo
import itertools

Hornwort = ['Anthoceros-agrestis-Bonn', 'Anthoceros-fusiformis', 'Anthoceros-punctatus', 'Anthoceros-agrestis-Oxford', 'Megaceros-flagellaris', 'Phaeomegaceros-chiloensis', 'Phymatoceros-phymatodes', 'Paraphymatoceros-pearsonii', 'Phaeoceros-carolinianus', 'Phaeoceros-sp', 'Phaeoceros-laevis-1032', 'Phaeoceros-laevis-902', 'Notothylas-orbicularis', 'Notothylas-ynanensis', 'Leiosporoceros-dussii']

Moss = ['Aerobryopsis-subdivergens', 'Cratoneuron-filicinum', 'Hygroamblystegium-varium', 'Hygrohypnum-luridum', 'Gollania-philippinensis', 'Hypnum-curvifolium', 'Calohypnum-plumiforme', 'Ectropothecium-obtusulum', 'Pleurozium-schreberi', 'Entodon-concinnus', 'Anomodon-attenuatus', 'Climacium-americanum', 'Rhytidiadelphus-subpinnatus', 'Brachythecium-laetum', 'Bryhnia-novae-angliae', 'Kindbergia-praelongum', 'Bryoandersonia-illecebra', 'Homaliodendron-scalpellifolium', 'Thamnobryum-sandei', 'Entodon-seductrix', 'Regmatodon-serrulatus', 'Calliergon-cordifolium', 'Catagonium-nitens-600', 'Catagonium-nitens-618', 'Catagonium-nitens-617', 'Lepyrodon-lagurus', 'Distichophyllum-collenchymatosum', 'Fontinalis-antipyretica', 'Fontinalis-sullivantii', 'Ptychomnion-cygnisetum', 'Hypopterygium-elatum', 'Hypopterygium-flavolimbatum', 'Orthotrichum-anomalum', 'Ulota-hutchinsiae', 'Racopilum-cuspidigerum', 'Mnium-hornum', 'Plagiomnium-ciliare', 'Pseudobryum-cinclidioides', 'Pohlia-nutans', 'Ptychostomum-knowltonii', 'Aulacomnium-androgynum', 'Aulacomnium-palustre', 'Bartramia-ithyphylla', 'Bartramia-mossmaniana', 'Philonotis-turneriana', 'Leptotheca-gaudichaudii', 'Pyrrhobryum-spiniforme', 'Splachnum-ampullaceum', 'Tayloria-subglabra', 'Aulacomnium-turgidum', 'Hedwigia-ciliata-1009', 'Hedwigia-ciliata-329', 'Archidium-alternifolium', 'Leucobryum-albidum', 'Leucobryum-bowringii', 'Grimmia-obtusifolia', 'Racomitrium-nitidulum', 'Niphotrichum-japonicum', 'Schistidium-cupulare', 'Ptychomitrium-wilsonii', 'Chorisodontium-acidophyllum', 'Dicranoloma-chilense', 'Dicranum-fulvum', 'Paraleucobryum-enerve', 'Barbula-amplexifolia', 'Gymnostomum-aurantiacum', 'Syntrichia-ruralis', 'Ceratodon-purpureus', 'Fissidens-javanicus', 'Schistostega-pennata', 'Pleuridium-subulatum', 'Bryoxiphium-norvegicum', 'Timmia-megapolitana', 'Encalypta-ciliata-1023', 'Encalypta-ciliata-1125', 'Gigaspermum-repens', 'Funaria-hygrometrica', 'Physcomitrella-patens', 'Physcomitrellopsis-africana', 'Diphyscium-fulvifolium', 'Buxbaumia-aphylla', 'Atrichum-angustatum', 'Pogonatum-microstomum', 'Pogonatum-subfuscatum', 'Polytrichastrum-ohioense', 'Polytrichum-strictum', 'Polytrichum-commune', 'Polytrichadelphus-magellanicus', 'Polytrichastrum-alpinum', 'Tetraphis-pellucida-1011', 'Tetraphis-pellucida-328', 'Andreaea-rupestris', 'Andreaea-wilsonii', 'Andreaeobryum-macrosporum', 'Takakia-lepidozioides', 'Sphagnum-fallax', 'Sphagnum-magellanicum', 'Sphagnum-girgensohnii', 'Sphagnum-palustre']

Liverwort = ['Acrobolbus-urvilleanus', 'Jungermannia-erectum', 'Mesoptychia-sp', 'Odontoschisma-sphagni', 'Plicanthus-hirtellus', 'Tetralophozia-filiformis', 'Scapania-nepalensis', 'Tritomaria-exsectiformis', 'Bazzania-tridens', 'Lepidozia-reptans', 'Herbertus-kurzii', 'Plagiochila-semidecurrens', 'Vetaforma-dusenii', 'Porella-caespitans-var-nipponica', 'Porella-chinensis', 'Porella-platyphylla', 'Ptilidium-pulcherrimum', 'Gackstroemia-magellanica', 'Acrolejeunea-sandvicensis', 'Frullania-moniliata', 'Metzgeria-furcata', 'Metzgeria-hamata', 'Pleurozia-purpurea', 'Fossombronia-cristula', 'Noteroclada-confluens', 'Pallavicinia-ambigua', 'Blasia-pusilla-349', 'Blasia-pusilla-890', 'Conocephalum-conicum', 'Cyathodium-cavernarum', 'Ricciocarpos-natans', 'Plagiochasma-appendiculatum', 'Riccia-sorocarpa', 'Marchantia-paleacea', 'Marchantia-quadrata', 'Marchantia-polymorpha', 'Sphaerocarpos-donnellii', 'Lunularia-cruciata', 'Haplomitrium-mnioides', 'Treubia-lacunosa']

Tracheophyte = ['Lycopodium-clavatum', 'Huperzia-asiatica', 'Diphasiastrum-complanatum', 'Selaginella-lepidophylla', 'Selaginella-kraussiana', 'Selaginella-moellendorffii', 'Selaginella-tamariscina', 'Isoetes-sinensis', 'Isoetes-taiwanensis', 'Dipteris-shenzhenensis', 'Alsophila-spinulosa', 'Adiantum-capillus-veneris', 'Ceratopteris-richardii', 'Azolla-filiculoides', 'Salvinia-cucullata', 'Marsilea-vestita', 'Torreya-grandis', 'Taxus-wallichiana', 'Taxus-chinensis', 'Sequoiadendron-giganteum', 'Ginkgo-biloba', 'Pinus-densiflora', 'Picea-mariana', 'Cycas-panzhihuaensis', 'Gnetum-montanum', 'Welwitschia-mirabilis', 'Amborella-trichopoda', 'Brasenia-schreberi', 'Nymphaea-colorata', 'Euryale-ferox', 'Ceratophyllum-demersum', 'Aristolochia-fimbriata', 'Piper-nigrum', 'Saururus-chinensis', 'Chloranthus-spicatus', 'Chloranthus-sessilifolius', 'Warburgia-ugandensis', 'Annona-cherimola', 'Magnolia-sieboldii', 'Liriodendron-chinense', 'Persea-americana', 'Chimonanthus-salicifolius', 'Acorus-gramineus', 'Wolffia-australiana', 'Acanthochlamys-bracteata', 'Dioscorea-alata', 'Angraecum-sesquipedale', 'Gloriosa-superba', 'Lilium-sargentiae', 'Oryza-sativa', 'Vriesea-erythrodactylon', 'Elaeis-guineensis', 'Pontederia-cordata', 'Musa-schizocarpa', 'Stephania-japonica', 'Protea-cynaroides', 'Tetracentron-sinense', 'Buxus-austroyunnanensis', 'Cornus-wilsoniana', 'Rhododendron-liliiflorum', 'Apocynum-pictum', 'Dracocephalum-rupestre', 'Ehretia-macrophylla', 'Eucommia-ulmoides', 'Solanum-lycopersicum', 'Nicotiana-benthamiana', 'Helwingia-omeiensis', 'Lactuca-sativa', 'Centella-asiatica', 'Triplostegia-glandulifera', 'Escallonia-herrerae', 'Malania-oleifera', 'Beta-vulgaris', 'Liquidambar-styraciflua', 'Vitis-rotundifolia', 'Glycine-max', 'Malus-doumeri', 'Quercus-variabilis', 'Citrullus-ecirrhosus', 'Tetraena-mongolica', 'Lagerstroemia-speciosa', 'Euscaphis-japonica', 'Populus-tremula', 'Averrhoa-carambola', 'Tripterygium-wilfordii', 'Gossypium-longicalyx', 'Arabidopsis-thaliana', 'Acer-yangbiense', 'Iodes-seguinii']

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    :This script is to curate representative four-taxon tree from a given gene tree
    :Usage python buildft.py buildft treefile -o treefile.ft
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('treefile', type=click.Path(exists=True))
@click.option('--output', '-o', default=None, show_default=True, help="file name of output")
def buildft(treefile,output):
    """
    :The tip labels of the input treefile need to be renamed as species names such as Anthoceros-agrestis-Bonn
    """
    Tree = Phylo.read(treefile,'newick')
    All_Clade_names = [clade.name for clade in Tree.find_clades() if clade.name]
    dis = {}
    clade_dict = {clade.name: clade for clade in Tree.find_clades() if clade.name}
    Hornwort_ = [c for c in Hornwort if c in clade_dict]
    Moss_ = [c for c in Moss if c in clade_dict]
    Liverwort_ = [c for c in Liverwort if c in clade_dict]
    Tracheophyte_ = [c for c in Tracheophyte if c in clade_dict]
    if len(Hornwort_) == 0:
        logging.info("No hornwort in the tree! Exit!")
        exit(0)
    if len(Moss_) == 0:
        logging.info("No moss in the tree! Exit!")
        exit(0)
    if len(Liverwort_) == 0:
        logging.info("No liverwort in the tree! Exit!")
        exit(0)
    if len(Tracheophyte_) == 0:
        logging.info("No tracheophyte in the tree! Exit!")
        exit(0)
    dist_cache = {}
    spname_clade_map = {}
    for i in Hornwort_: spname_clade_map.update({i:"Hornwort"})
    for i in Moss_: spname_clade_map.update({i:"Moss"})
    for i in Liverwort_: spname_clade_map.update({i:"Liverwort"})
    for i in Tracheophyte_: spname_clade_map.update({i:"Tracheophyte"})
    for c1, c2 in itertools.combinations(clade_dict.values(), 2):
        mrca = Tree.common_ancestor(c1, c2)
        d = abs(mrca.distance(c1) - mrca.distance(c2))
        dist_cache[(c1.name, c2.name)] = d
        dist_cache[(c2.name, c1.name)] = d
    for i in Hornwort_:
        for j in Moss_:
            for k in Liverwort_:
                for l in Tracheophyte_:
                    names = [i, j, k, l]
                    clades_id = ", ".join(names)
                    if any(dist_cache[(c1,c2)]<=0 for c1,c2 in itertools.combinations(names,2)):
                        continue
                    dis_tmp = sum(dist_cache[(c1,c2)] for c1,c2 in itertools.combinations(names,2))
                    dis[clades_id] = dis_tmp
    best_4_sp = sorted(dis.items(),key=lambda x:x[1])[0][0]
    best_4_sp = best_4_sp.split(', ')
    for tip in Tree.get_terminals():
        if tip.name not in best_4_sp:
            Tree.prune(tip)
    for tip in Tree.get_terminals():
        tip.name = spname_clade_map[tip.name]
    Phylo.write(Tree,output,format='newick')

if __name__ == '__main__':
        cli()
