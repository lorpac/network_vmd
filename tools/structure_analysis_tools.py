import os
from Bio.PDB import PDBList
import networkx as nx
import biographs as bg
import pandas as pd
import numpy as np
import pickle
import tools.amino_acids_conversion as aaconv
import glob
from biopandas.pdb import PandasPdb


def download_pdb(pdb_id, pdbs_path):
    """Downloads a pdb file
    
    Parameters
    ----------
    pdb_id : str
        pdb id
    pdbs_path: str, optional
        path of folder containing the pdbs (default is "pdbs")
    """
    pdbl = PDBList(obsolete_pdb=True)
    pdbl.download_pdb_files(pdb_id, file_format="pdb", pdir=pdbs_path)


def get_pdb_path(pdb_id, pdbs_path="pdbs"):
    """Gets the location of a pdb file
    
    Parameters
    ----------
    pdb_id : str
        pdb id
    pdbs_path: str, optional
        path of folder containing the pdbs (default is "pdbs")
    
    Returns
    -------
    str
        pdb file path
    str
        True if the pdb has been downloaded
    """

    downloaded = False

    if not pdbs_path:
        pdbs_path = "pdbs"

    abs_path = os.path.join("data", pdbs_path)
    abs_file_path = os.path.join(abs_path, pdb_id + ".*")

    if len(glob.glob(abs_file_path)) == 0:
        abs_file_path = os.path.join(abs_path, "pdb" + pdb_id + ".*")
        downloaded = True

        if len(glob.glob(abs_file_path)) == 0:
            os.makedirs(abs_path, exist_ok=True)
            download_pdb([pdb_id], abs_path)

        else:
            pdb_id = "pdb" + pdb_id
            abs_file_path = os.path.join(abs_path, pdb_id + ".*")
    pdb_path = glob.glob(abs_file_path)[0]

    return pdb_path, downloaded


def remove_hydrogens(pdb_file):
    """Remove hydrogens from a pdb file
    
    Parameters
    ----------
    pdb_file : str
        pdb file path
    
    Returns
    -------
    list
        old pdb file lines
    list
        new odb file lines
    """

    with open(pdb_file, "r") as f:
        pdb = [line.rsplit("\n")[0] for line in f]
    pdb_new = []
    for line in pdb:
        if line[0:4] == "ATOM" and line[13] == "H":
            pass
        else:
            pdb_new.append(line)
    with open(pdb_file, "w") as f:
        f.writelines([line + "\n" for line in pdb_new])

    return pdb, pdb_new


def list_relations(dim="all"):
    """Returns the list of relations to analyze
    
    Parameters
    ----------
    dim : str, optional
        "1D", "2D", "1-2D", "3-4D", "all" or "" (default is "all")
        "" is interpreted as "all"
    
    Returns
    -------
    list
        list of relations to analyze
    """

    if dim == "2D":
        rel_list = ["2D2", "2D3", "2D4"]
    elif dim == "1D":
        rel_list = ["1D"]
    elif dim == "1-2D":
        rel_list = ["1D", "2D2", "2D3", "2D4"]
    elif dim == "3-4D":
        rel_list = ["3D", "4D"]
    elif dim == "all" or dim == "" :
        rel_list = ["1D", "2D2", "2D3", "2D4", "3D", "4D"]
    else:
        rel_list = [dim]

    return rel_list


def assign_secondary_structure(pdb):
    """Returns the secondary structure elements of a pdb
    
    Parameters
    ----------
    pdb : str
        pdb file path
    
    Returns
    -------
    dict
        dictionary of secondary structure elements
    """

    ppdb = PandasPdb().read_pdb(pdb)

    secondary_structure = {}

    helices_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "HELIX"][
        "entry"
    ]
    for helix in helices_from_pdb:
        identifier_h = helix[5:8].strip()
        initial_chain_h = helix[13].strip()
        initial_pos_h = helix[16:19].strip()
        final_pos_h = helix[28:31].strip()
        for i in range(int(initial_pos_h), int(final_pos_h) + 1):
            secondary_structure[initial_chain_h + str(i)] = (
                "helix" + identifier_h + "-" + initial_chain_h
            )

    sheets_from_pdb = ppdb.df["OTHERS"][ppdb.df["OTHERS"]["record_name"] == "SHEET"][
        "entry"
    ]
    for sheet in sheets_from_pdb:
        identifier_s = sheet[6:8].strip()
        initial_chain_s = sheet[15].strip()
        initial_pos_s = sheet[17:20].strip()
        final_pos_s = sheet[28:31].strip()
        for i in range(int(initial_pos_s), int(final_pos_s) + 1):
            secondary_structure[initial_chain_s + str(i)] = (
                "sheet" + identifier_s + "-" + initial_chain_s
            )

    mol = bg.Pmolecule(pdb)
    net = mol.network()

    residues_type = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname
        res_pos = residue.parent.id + str(residue.id[1])
        residues_type[res_pos] = res_type

    residues = list(net.nodes)  # assume they are ordered
    last_structure = None
    last_chain = None
    i = 0
    for residue in residues:
        chain = residue[0]
        try:
            structure = secondary_structure[residue]
            if structure != last_structure:
                i += 1
        except KeyError:
            if chain != last_chain:
                i += 1
            structure = "loop" + str(i)
            secondary_structure[residue] = structure
        last_structure = structure
        last_chain = chain

    return secondary_structure


def get_neighbor_structure_relation(secondary_structure, u, v):
    """Returns the relation (1D, 2D, 3D, 4D) between to neighboring nodes
    
    Parameters
    ----------
    secondary_structure : dict
        dictionary of secondary structure elements
    u: str
        node label
    v: str
        node label
    
    Returns
    -------
    str
        relation (1D, 2D, 3D, 4D)
    """

    chain_u = u[0]
    chain_v = v[0]
    pos_u = u[1::]
    pos_v = v[1::]
    struct_u = secondary_structure[u]
    struct_v = secondary_structure[v]

    if chain_u == chain_v:
        dist = np.abs(int(pos_u) - int(pos_v))
        if dist == 1:
            relation = "1D"
        elif struct_u == struct_v:
            if dist < 5:
                relation = "2D" + str(dist)
            else:
                relation = "3D"
        else:
            relation = "3D"
    else:
        relation = "4D"

    return relation


def create_aa_network(
    pdb_id,
    rel_list,
    selected_positions=None,
    cutoff=5,
    pdbs_path="pdbs",
):
    """Creates the amino acid network from a pdb id.
        
    Parameters
    ----------
    pdb_id : str
        pdb id of the protein
    rel_list: list
        list of relation (1D, 2D, 3D, 4D) to consider.
    folder_path: str
        path of the output folder
    selected_positions: None or list, optional
        list of sequence positions to consider. If None, all positions are considered (default is None)
    cutoff: int, optional
        cutoff threshold for the connection of nodes in the amino acid network (dafault is 5).
    pdbs_path: str, optional
        path of the pdb files folder (default is "pdbs")
    Returns
    -------
    Graph
        NetworkX Graph object
    dict
        labels dictionary
    """

    pdb, downloaded = get_pdb_path(pdb_id, pdbs_path)

    # initialize database to report
    database = []

    mol = bg.Pmolecule(pdb)
    net = mol.network(cutoff=cutoff)

    # take only selected positions:
    if selected_positions:
        for node in list(net.nodes):
            pos = int(node[1::])
            if pos not in selected_positions:
                net.remove_node(node)
    else:
        positions = [int(node[1::]) for node in list(net.nodes)]
        pos_min = min(positions)
        pos_max = max(positions)
        selected_positions = range(pos_min, pos_max + 1)

    secondary_structure = assign_secondary_structure(pdb)

    labels = {}
    residues_dict = {}
    for residue in mol.model.get_residues():
        res_type = residue.resname.strip()
        if len(res_type) > 1:
            res_type = aaconv.three2one(res_type)
        res_chain = residue.parent.id
        res_pos = str(residue.id[1])
        labels[res_chain + res_pos] = f"{res_type}{res_pos}:{res_chain}"
        residues_dict[res_chain + res_pos] = res_type
    
    for residue in mol.model.get_residues():
        node_name = residue.parent.id + str(residue.id[1])
        deg = nx.degree(net, node_name)
        if deg == 0:
            net.remove_node(node_name)
            _ = labels.pop(node_name)
        else:
            seqpos = residue.id[1]
            if seqpos not in selected_positions:
                net.remove_node(node_name)
                _ = labels.pop(node_name)
            structure = secondary_structure[node_name]

            for neighbor in list(nx.neighbors(net, node_name)):
                edge_weight = nx.edges(net)[(node_name, neighbor)]["weight"]
                relation = get_neighbor_structure_relation(
                    secondary_structure, node_name, neighbor
                )
                # select only edges of desired relations
                if relation not in rel_list:
                    net.remove_edge(neighbor, node_name)

            # check if the residue became of degree zero:
            deg = nx.degree(net, residue.parent.id + str(residue.id[1]))

            if deg == 0:
                net.remove_node(node_name)
                _ = labels.pop(node_name)

    nx.relabel_nodes(net, labels, copy=False)

    return net