# Visualize amino acid networks in VMD

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
        --pdbs_path TEXT           path of folder containing the pdbs within the
                                    data folder.  [default: pdbs]
        --help                     Show this message and exit.
        ```

        Example:
        ```console
        C:\Users\Lorenza\Documents\VMDscript>py create_net.py 1f41
        ```
    - **Perturbation Network:** Run the `create_perturbation_net.py` script:

        ```console
        Usage: create_perturbation_net.py [OPTIONS] PDB_ID PDB_ID_REF

        Options:
        --threshold INTEGER        difference threshold do connect residues in the
                                    perturbation network  [default: 0]
        --select_positions TEXT    if True, create a list of sequence positions to
                                    consider (from --start and --stop options)
                                    [default: False]
        --range_positions TEXT...  if not None and --select_positions is True, first
                                    and last positions of sequence to consider. If
                                    None or --select_positions is False, consider the
                                    whole sequence.
        --dim TEXT                 dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or
                                    all).  [default: all]
        --cutoff INTEGER           cutoff threshold for the connection of nodes in
                                    the amino acid network.  [default: 5]
        --output TEXT              output file name. If None, the output is
                                    'netdimpdb_id.p'.
        --output_folder TEXT       output folder. Default is 'networks'.
        --pdbs_path TEXT           path of folder containing the pdbs within the
                                    data folder (default is 'pdbs')
        --help                     Show this message and exit.
        ```

        Example:
        ```console
        C:\Users\Lorenza\Documents\VMDscript> py .\create_perturbation_net.py 3djz 1f41 --select_positions True --range_positions 11 124  --threshold 4
        ```

2. Run the `network_vmd.py` script:

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
    -nc NC                the network has no edge color attribute
    -ntodraw NTODRAW [NTODRAW ...]
                            draw only a list of nodes
    -norm NORM            changes the normalizaton factor
    ```

    Examples:
    ```console
    C:\Users\Lorenza\Documents\VMDscript> py network_vmd.py data/pdbs/pdb1f41.ent networks/net1f41.p tcl_files/net1f41.tcl -nc 1
    ```

    ```console
    C:\Users\Lorenza\Documents\VMDscript> py network_vmd.py data/pdbs/pdb3djz.ent networks/pert_net3djz.p tcl_files/pert_net3djz.tcl
    ```

    NOTE: create_net.py has saved a copy of the pdb file in data/pdbs/

3. Open the pdb file with VMD.

4. In the Tk console, navigate to your folder and run source output.tcl .

    Example:
    ```console
    Main console display active (Tcl8.5.6 / Tk8.5.6)
    (VMD) 1 % cd ../../../Users/Lorenza/Documents/VMDscript/
    >Main< (VMDscript) 2 % source tcl_files/net1f41.tcl
    ```