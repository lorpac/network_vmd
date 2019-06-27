import networkx as nx
import structure_analysis_tools as sa
import amino_acids_conversion as aac
import os

import click


@click.command()
@click.argument('pdb_id')
@click.option('--dim', default='all', help="dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or all). Default is 'all'. ")
@click.option('--output', default=None, help="output file name. If None (default), the output is 'netdimpdb_id.p'.")
@click.option('--output_folder', default="networks", help="output folder. Default is 'networks'.")
def create_net(pdb_id, dim, output, output_folder):
    if dim == "all": dim = ""
    rel_list = sa.list_relations(dim)
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

    if not output:
        output = dim + "net" + pdb_id + ".p"

    os.makedirs(output_folder, exist_ok=True)

    nx.write_gpickle(net, os.path.join(output_folder, output))

if __name__ == '__main__':
    create_net()