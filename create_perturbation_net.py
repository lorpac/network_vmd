import networkx as nx
import tools.structure_analysis_tools as sa
import os
import click


@click.command()
@click.argument("pdb_id")
@click.argument("pdb_id_ref")
@click.option(
    "--threshold",
    default=0,
    help="difference threshold do connect residues in the perturbation network",
    show_default=True
)
@click.option(
    "--select_positions",
    default=False,
    help="if True, create a list of sequence positions to consider (from --start and --stop options)",
    show_default=True
)
@click.option(
    "--range_positions",
    default=None,
    nargs=2,
    help="if not None and --select_positions is True, first and last positions of sequence to consider. If None or --select_positions is False, consider the whole sequence.",
    show_default=True
)
@click.option(
    "--dim",
    default="all",
    help="dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or all).",
    show_default=True
)
@click.option(
    "--cutoff",
    default=5,
    help="cutoff threshold for the connection of nodes in the amino acid network.",
    show_default=True
)
@click.option(
    "--output",
    default=None,
    help="output file name. If None, the output is 'netdimpdb_id.p'.",
    show_default=True
)
@click.option(
    "--output_folder", default="networks", help="output folder. Default is 'networks'.",
    show_default=True
)
@click.option(
    "--pdbs_path",
    default="pdbs",
    help="path of folder containing the pdbs.",
    show_default=True
)
def create_net(
    pdb_id,
    pdb_id_ref,
    threshold,
    select_positions,
    range_positions,
    dim,
    cutoff,
    output,
    output_folder,
    pdbs_path,
):

    if dim == "all":
        dim = ""

    rel_list = sa.list_relations(dim)

    if select_positions:
        start, stop = (int(x) for x in range_positions)
        selected_positions = list(range(start, stop + 1))
    else:
        selected_positions = None
        

    net, labels = sa.create_perturbation_network(
        pdb_id,
        pdb_id_ref,
        threshold,
        rel_list,
        selected_positions=selected_positions,
        cutoff=cutoff,
        pdbs_path=pdbs_path,
    )

    nx.relabel_nodes(net, labels, copy=False)

    if not output:
        output = dim + "pert_net" + pdb_id + ".p"

    os.makedirs(output_folder, exist_ok=True)

    nx.write_gpickle(net, os.path.join(output_folder, output))


if __name__ == "__main__":
    create_net()
