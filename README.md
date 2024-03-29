# Visualize amino acid networks in VMD

The following Python packages are required:
- [biographs](https://github.com/rodogi/biographs)
- [Biopython](https://biopython.org/wiki/Download)
- [NetworkX](https://networkx.github.io/)
- [os](https://docs.python.org/3/library/os.html)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://www.numpy.org/)
- [pickle](https://docs.python.org/3/library/pickle.html)
- [click](https://click.palletsprojects.com/en/7.x/)

The results can be visualized with the [VMD software](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)

## Instructions:

1. - **Amino Acid network:** Run the `create_net.py` script:

        ```console
        Usage: create_net.py [OPTIONS] PDB_ID

        Options:
        --select_positions TEXT    if True, create a list of sequence positions to
                                    consider (from --start and --stop options)
                                    [default: False]
        --range_positions TEXT...  if not None and --select_positions is True, first
                                    and last positions of sequence to consider. If
                                    None or --select_positions is False, consider the
                                    whole sequence.
        --dim TEXT                 dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or
                                    all)  [default: all]
        --cutoff INTEGER           cutoff threshold for the connection of nodes in
                                    the amino acid network  [default: 5]
        --output TEXT              output file name. If None (default), the output
                                    is 'netdimpdb_id.p'.
        --output_folder TEXT       output folder.  [default: networks]
        --pdbs_path TEXT           path of folder containing the pdbs.  [default:
                                    pdbs]
        --help                     Show this message and exit.
        ```

        Example:
        ```console
        C:\Users\Lorenza\Documents\network_to_vmd>py create_net.py 1f41
        ```
    - **Perturbation Network:** Run the `create_perturbation_net.py` script:

        ```console
        Usage: create_perturbation_net.py [OPTIONS] PDB_ID PDB_ID_REF

        Options:
        --threshold INTEGER         difference threshold do connect residues in the
                                    perturbation network  [default: 0]
        --select_positions TEXT     if True, create a list of sequence positions to
                                    consider (from --start and --stop options)
                                    [default: False]
        --range_positions TEXT...   if not None and --select_positions is True,
                                    first and last positions of sequence to
                                    consider. If None or --select_positions is
                                    False, consider the whole sequence.
        --dim TEXT                  dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or
                                    all).  [default: all]
        --cutoff INTEGER            cutoff threshold for the connection of nodes in
                                    the amino acid network.  [default: 5]
        --output TEXT               output file name. If None, the output is
                                    'netdimpdb_id.p'.
        --output_folder TEXT        output folder. Default is 'networks'.  [default:
                                    networks]
        --pdbs_path TEXT            path of folder containing the pdbs.  [default:
                                    pdbs]
        --color_negative_edge TEXT  color of negative edges (wij_ref > wij).
                                    [default: red]
        --color_positive_edge TEXT  color of positive edges (wij > wij_ref).
                                    [default: green]
        --help                      Show this message and exit.
        ```

        Example:
        ```console
        C:\Users\Lorenza\Documents\network_to_vmd> py .\create_perturbation_net.py 3djz 1f41 --select_positions True --range_positions 11 124  --threshold 4
        ```

2. Run the `network_vmd.py` script (adapted from Aria Gheeraert's [PerturbationNetworkAnalysis](https://github.com/agheeraert/PerturbationNetworkAnalysis)):

    ```console
    usage: network_vmd.py [-h] [-nc NC] [-ntodraw NTODRAW [NTODRAW ...]]
                      [-norm NORM]
                      pdb_path net_path output

    file paths

    positional arguments:
    pdb_path              pdb file
    net_path              network file
    output                output file

    optional arguments:
    -h, --help            show this help message and exit
    -nc 0 or 1                the network has no edge color attribute
    -ntodraw NTODRAW [NTODRAW ...]
                            draw only a list of nodes
    -norm NORM            changes the normalizaton factor
    -weighted            use edge weights for thickness
    -color                     color (ignored if -nc = False, default: red)
    ```

    Examples:
    ```console
    C:\Users\Lorenza\Documents\network_to_vmd> py network_vmd.py pdbs/pdb1f41.ent networks/net1f41.p tcl_files/net1f41.tcl -nc 1
    ```

    ```console
    C:\Users\Lorenza\Documents\network_to_vmd> py network_vmd.py pdbs/pdb3djz.ent networks/pert_net3djz_4.p tcl_files/pert_net3djz_4.tcl
    ```

    NOTE: create_net.py has saved a copy of the pdb file in pdbs/

3. Open the pdb file with VMD.

4. In the Tk console of VMD, navigate to the project folder and run `source output.tcl` .

    Example:
    ```console
    Main console display active (Tcl8.5.6 / Tk8.5.6)
    (VMD) 1 % cd ../../../Users/Lorenza/Documents/network_to_vmd/
    >Main< (network_to_vmd) 2 % source tcl_files/net1f41.tcl
    ```
