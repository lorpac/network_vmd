# Visualize amino acid networks in VMD

1. Run the `create_net.py` script:

    ```console
    Usage: create_net.py [OPTIONS] PDB_ID

    Options:
    --selected_positions TEXT  list of sequence positions to consider. If None,
                                all positions are considered (default is None)
    --dim TEXT                 dimensions to consider (1D, 2D, 3D, 4D, 3-4D, or
                                all). Default is 'all'.
    --cutoff INTEGER           cutoff threshold for the connection of nodes in
                                the amino acid network (dafault is 5)
    --output TEXT              output file name. If None (default), the output
                                is 'netdimpdb_id.p'.
    --output_folder TEXT       output folder. Default is 'networks'.
    --help                     Show this message and exit.
    ```

    Example:
    ```console
        C:\Users\Lorenza\Documents\VMDscript>py create_net.py 1f41
    ```

2. Run the `network_vmd.py` script:
then follow Aria's instructions (see below):

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

    Example:
    ```console
        C:\Users\Lorenza\Documents\VMDscript> py network_vmd.py data/pdbs/pdb1f41.ent networks/net1f41.p tcl_files/net1f41.tcl -nc 1
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