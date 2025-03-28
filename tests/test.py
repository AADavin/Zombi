import os
import filecmp
from typing import Dict



def fill_template(template_path: str, output_path: str, replacements: Dict[str, str]):
    with open(template_path, 'r') as file:
        content = file.read()
    for key, value in replacements.items():
        if key == 'ASSORTATIVE_TRANSFER_FIELD':
            # Special handling for boolean replacement
            content = content.replace(key, 'True')
        else:
            content = content.replace(key, str(value))
    with open(output_path, 'w') as file:
        file.write(content)

def run_zombi(mode: str, param_file: str, output_folder: str):
    command = f"python3 Zombi.py {mode} {param_file} {output_folder}"
    os.system(command)

def compare_folders(folder1: str, folder2: str) -> bool:
    comparison = filecmp.dircmp(folder1, folder2)
    if comparison.diff_files:
        print(f"Files that differ: {comparison.diff_files}")
        return False
    return True


def test_reproducibility():
    # Set up paths
    species_template = './Parameters/SpeciesTreeParametersTemplate.tsv'
    genome_template = './Parameters/GenomeParametersTemplate.tsv'
    os.makedirs("./Test_Outputs", exist_ok=True)
    output_folder_1 = './Test_Outputs/Output_folder_1'
    output_folder_2 = './Test_Outputs/Output_folder_2'
    test_parameters_dir = './Parameters/TestParameters'
    if os.path.exists(test_parameters_dir):
        # remove the directory if it exists
        os.system(f'rm -rf {test_parameters_dir}')
    # Create test parameters directory
    os.makedirs(test_parameters_dir, exist_ok=True)

    replacements = {
        'DUPLICATION_FIELD': '0.00',
        'TRANSFER_FIELD': '1',
        'LOSS_FIELD': '0.00',
        'ASSORTATIVE_TRANSFER_FIELD': 'True',
        'ALPHA_FIELD': '5'
    }

    # Fill templates
    fill_template(species_template, f'{test_parameters_dir}/SpeciesTreeParameters.tsv', replacements)
    fill_template(genome_template, f'{test_parameters_dir}/GenomeParameters.tsv', replacements)

    # Run Zombi twice
    os.makedirs(output_folder_1, exist_ok=True)
    os.makedirs(output_folder_2, exist_ok=True)
    run_zombi('T', f'{test_parameters_dir}/SpeciesTreeParameters.tsv', output_folder_1)
    run_zombi('Gm', f'{test_parameters_dir}/GenomeParameters.tsv', output_folder_1)

    run_zombi('T', f'{test_parameters_dir}/SpeciesTreeParameters.tsv', output_folder_2)
    run_zombi('Gm', f'{test_parameters_dir}/GenomeParameters.tsv', output_folder_2)

    # Compare outputs
    if compare_folders(output_folder_1, output_folder_2):
        print('The outputs are identical! ✅')
    else:
        print('The outputs are different! ❌')


test_reproducibility()