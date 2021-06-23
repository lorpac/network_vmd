# adapted from Aria Gheeraert's PerturbationNetworkAnalysis
# (https://github.com/agheeraert/PerturbationNetworkAnalysis)

import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import argparse

parser = argparse.ArgumentParser(description='file paths')
parser.add_argument('pdb_path',  type=str,
                    help='pdb file')
parser.add_argument('net_path',  type=str,
                    help='network file')
parser.add_argument('output',  type=str,
                    help='output file')
parser.add_argument('-nc',  type=bool, default=False,
                    help='the network has no edge color attribute')
parser.add_argument('-weighted',  action='store_true', default=False,
                    help='use edge weights')
parser.add_argument('-ntodraw', type=str, nargs='+', default=None,
                    help="draw only a list of nodes")
parser.add_argument('-neighborhood', type=str, default=None,
                    help="draw neighborhood of the node")
parser.add_argument('-norm', type=float, default=1.5,
                    help="changes the normalizaton factor")
parser.add_argument('-color', type=str, default='red',
                    help="color (ignored if nc = True)")

args = parser.parse_args()

if args.ntodraw:
    if 'range' in args.ntodraw:
        ntodraw = [str(i) for i in range(int(args.ntodraw[1]), int(args.ntodraw[2]))]
    else:
        ntodraw = args.ntodraw
else:
    ntodraw = None
three2one = dict(zip(aa3, aa1))
three2one['5CS'] = 'C'
A = nx.read_gpickle(args.net_path)
structure = PDBParser().get_structure('X', args.pdb_path)[0]
color = not args.nc
div = max(nx.get_edge_attributes(A, 'weight').values())/args.norm

node2CA = {}
for atom in structure.get_atoms():
    if atom.id == 'CA':
        residue = atom.parent
        if residue.resname in three2one:
                coords = str(atom.coord[0]) + ' ' + str(atom.coord[1]) + ' ' + str(atom.coord[2])
                node2CA[three2one[residue.resname]+str(residue.id[1])+':'+residue.parent.id] = coords

with open(args.output, 'w') as output:
        output.write('draw delete all \n')
        if color:
            output.write('draw color green \n')
        else:
            output.write('draw color %s \n' %args.color)
        previous = None
        for u, v in A.edges():
            if color:
                c = A.get_edge_data(u, v)['color']
                if previous != c:
                    output.write(f'draw color {c} \n')
                    previous = c
                if ntodraw:
                    if u[1:-2] in ntodraw and v[1:-2] in ntodraw:
                        if args.weighted:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                        else:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(0.5)+' \n')
                else:
                    if args.weighted:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                    else:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(0.5)+' \n')
            else:
                if ntodraw:
                    if u[1:-2] in ntodraw and v[1:-2] in ntodraw:
                        if args.weighted:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                        else:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(0.5)+' \n')
                elif args.neighborhood:
                    if u == args.neighborhood or v == args.neighborhood:
                        if args.weighted:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                        else:
                            output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(0.5)+' \n')

                else:
                    if args.weighted:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(A.get_edge_data(u, v)['weight']/div)+' \n')
                    else:
                        output.write('draw cylinder { ' + str(node2CA[u]) + ' '+ ' } ' + '{ ' + str(node2CA[v]) + ' '+ ' } radius '+str(0.5)+' \n')

        output.write('draw color 2 \n')

        for u in A.nodes():
            if ntodraw:
                if u in ntodraw[1:-2]:
                    output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(args.norm)+' \n')
            elif args.neighborhood:
                if u == args.neighborhood or u in A.neighbors(u):
                    output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(args.norm)+' \n')
            else:
                output.write('draw sphere { ' + str(node2CA[u]) + ' '+ ' } radius '+str(args.norm)+' \n')
