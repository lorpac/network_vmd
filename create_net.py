import networkx as nx
import tools.structure_analysis_tools as sa
import os
import click


@click.command()
@click.argument("pdb_id")
@click.option(
    "--selected_positions",
    default=None,
    help="list of sequence positions to consider. If None, all positions are considered (default is None)",
)
@click.option(
    "--dim",
    default="all",
    help="dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or all). Default is 'all'. ",
)
@click.option(
    "--cutoff",
    default=5,
    help="cutoff threshold for the connection of nodes in the amino acid network (dafault is 5)",
)
@click.option(
    "--output",
    default=None,
    help="output file name. If None (default), the output is 'netdimpdb_id.p'.",
)
@click.option(
    "--output_folder", default="networks", help="output folder. Default is 'networks'."
)
@click.option(
    "--pdbs_path",
    default="pdbs",
    help="path of folder containing the pdbs within the data folder (default is 'pdbs')",
)
def create_net(
    pdb_id, selected_positions, dim, cutoff, output, output_folder, pdbs_path
):

    if dim == "all":
        dim = ""

    rel_list = sa.list_relations(dim)

    net = sa.create_aa_network(
        pdb_id,
        rel_list,
        selected_positions=selected_positions,
        cutoff=cutoff,
        pdbs_path=pdbs_path,
    )

    if not output:
        output = dim + "net" + pdb_id + ".p"

    os.makedirs(output_folder, exist_ok=True)

    nx.write_gpickle(net, os.path.join(output_folder, output))


if __name__ == "__main__":
    create_net()
