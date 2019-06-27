import networkx as nx
import sys
import structure_analysis_tools as sa
import amino_acids_conversion as aac

_, pdb_id = sys.argv
rel_list = sa.list_relations("all")

net, db_1, db_2, mol, downloaded = sa.create_aa_network(
    pdb_id,
    rel_list,
    folder_path=None,
    save_csv=False,
)

node_labels = {}

for n in net.nodes:
    info = db_2[db_2["Position"] == n]
    chain = n[0]
    pos = n[1::]
    typeaa = info["Type of residue"].item()
    typeaa = aac.three2one(typeaa)

    label = f"{typeaa}{pos}:{chain}"

    node_labels[n] = label

nx.relabel_nodes(net, node_labels, copy=False)

nx.write_gpickle(net, pdb_id + ".p")