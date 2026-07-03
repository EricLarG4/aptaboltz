from pymol import cmd
import os
import numpy as np
import pandas as pd


## --- Prediction processing ---

def process_experiment(project, model, experiment):
    print(f"Processing experiment: {project}/{model}/{experiment}")

    # file path definitions
    cif_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{cif_path}"

    # list cif files in the specified experiment directory
    cif_files = [f for f in os.listdir(full_path) if f.endswith(".cif")]

    # sort cif_files by filenames naturally to keep best models first
    cif_files.sort(key=lambda x: int(x.split("_")[-1].split(".")[0]))

    print(f"Found CIF files: {cif_files}")

    cmd.reinitialize()

    # loading file with PyMOL
    for cif_file in cif_files:
        print(f"Loading CIF file: {full_path + cif_file}")
        cmd.load(full_path + cif_file)

    # align, join states and save multi-state file
    cmd.alignto()
    cmd.join_states(f"{model}", "input_model*")
    cmd.delete("input_model*")
    cmd.save(f"{experiment_path}{model}_{experiment}.cif", state=0)

    print(f"Saved multi-state file: {experiment_path}{model}_{experiment}.cif")
    cmd.reinitialize()
    cmd.load(f"{experiment_path}{model}_{experiment}.cif")
    cmd.load("python/base.py")
    cmd.save(f"{experiment_path}{model}_{experiment}.pse")

## --- Confidence processing ---

def process_confidence(project, model, experiment):
    print(f"Processing confidence for experiment: {project}/{model}/{experiment}")

    # file path definitions
    json_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{json_path}"

    # list JSON files in the specified experiment directory
    json_files = [f for f in os.listdir(full_path) if f.endswith(".json")]

    print(f"Found JSON files: {json_files}")

    # process each JSON file and create a pandas dataframe
    dataframes = []
    for json_file in json_files:
        json_file_path = os.path.join(full_path, json_file)
        df = pd.read_json(json_file_path)

        # clean dict-like columns: extract the value from {'0': value} -> value
        for col in ['pair_chains_iptm', 'pair_chains_pae']:
            if col in df.columns:
                df[col] = df[col].apply(
                    lambda x: list(x.values())[0] if isinstance(x, dict) else x
                )

        dataframes.append(df)

    # concatenate all dataframes into a single dataframe
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(f"{experiment_path}{model}_{experiment}_confidence.csv", index=False)
        print(f"Saved combined CSV file: {experiment_path}{model}_{experiment}.csv")

def process_pae(project, model, experiment):
    print(f"Processing PAE for experiment: {project}/{model}/{experiment}")

    # file path definitions
    pae_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{pae_path}"

    # list PAE JSON files in the specified experiment directory
    pae_files = [f for f in os.listdir(full_path) if f.endswith(".npz") and f.startswith("pae_input_model_")]

    print(f"Found PAE files: {pae_files}")

    # process each PAE file and create a pandas dataframe
    dataframes = []
    for pae_file in pae_files:
        pae_file_path = os.path.join(full_path, pae_file)
        # Load the .npz file and extract the PAE matrix
        with np.load(pae_file_path) as data:
            pae_matrix = data['pae']
            df = pd.DataFrame(pae_matrix)
            df['model'] = pae_file.split('_')[-1].split('.')[0]  # Extract model number from filename
            dataframes.append(df)
    
    # concatenate all dataframes into a single dataframe
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(f"{experiment_path}{model}_{experiment}_pae.csv", index=False)
        print(f"Saved combined PAE CSV file: {experiment_path}{model}_{experiment}_pae.csv")

def process_pde(project, model, experiment):
    print(f"Processing PDE for experiment: {project}/{model}/{experiment}")

    # file path definitions
    pde_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{pde_path}"

    # list PDE JSON files in the specified experiment directory
    pde_files = [f for f in os.listdir(full_path) if f.endswith(".npz") and f.startswith("pde_input_model_")]

    print(f"Found PDE files: {pde_files}")

    # process each PDE file and create a pandas dataframe
    dataframes = []
    for pde_file in pde_files:
        pde_file_path = os.path.join(full_path, pde_file)
        # Load the .npz file and extract the PDE matrix
        with np.load(pde_file_path) as data:
            pde_matrix = data['pde']
            df = pd.DataFrame(pde_matrix)
            df['model'] = pde_file.split('_')[-1].split('.')[0]  # Extract model number from filename
            dataframes.append(df)
    
    # concatenate all dataframes into a single dataframe
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(f"{experiment_path}{model}_{experiment}_pde.csv", index=False)
        print(f"Saved combined PDE CSV file: {experiment_path}{model}_{experiment}_pde.csv")

def process_plddt(project, model, experiment):
    print(f"Processing pLDDT for experiment: {project}/{model}/{experiment}")

    # file path definitions
    plddt_path = "boltz_results_input/predictions/input/"
    experiment_path = f"{project}/{model}/{experiment}/"
    full_path = f"{experiment_path}{plddt_path}"

    # list pLDDT JSON files in the specified experiment directory
    plddt_files = [f for f in os.listdir(full_path) if f.endswith(".npz") and f.startswith("plddt_input_model_")]

    print(f"Found pLDDT files: {plddt_files}")

    # process each pLDDT file and create a pandas dataframe
    dataframes = []
    for plddt_file in plddt_files:
        plddt_file_path = os.path.join(full_path, plddt_file)
        # Load the .npz file and extract the pLDDT array
        with np.load(plddt_file_path) as data:
            plddt_array = data['plddt']
            df = pd.DataFrame(plddt_array, columns=['pLDDT'])
            df['model'] = plddt_file.split('_')[-1].split('.')[0]  # Extract model number from filename
            dataframes.append(df)
    
    # concatenate all dataframes into a single dataframe
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        combined_df.to_csv(f"{experiment_path}{model}_{experiment}_plddt.csv", index=False)
        print(f"Saved combined pLDDT CSV file: {experiment_path}{model}_{experiment}_plddt.csv")