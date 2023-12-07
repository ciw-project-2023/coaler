#!/usr/bin/env python3
import os
import subprocess
import time
from dataclasses import dataclass

import sys
import tempfile
from os import listdir
from rdkit import Chem
from rdkit.Chem import AllChem

target = 'Release'
verbose = True
num_conformers = 20
threads = 16
pwd = os.getcwd()
nbse_dir = 'nbse'
coaler_bin = 'build/{}/src/coaler'.format(target)
def flatten(l):
    return [item for sublist in l for item in sublist]

def sdf_to_smiles(file: str) -> list[Chem.Mol]:
    smiles = []
    suppl = Chem.SDMolSupplier(file)
    for mol in suppl:
        smiles.append(mol)

    return smiles


def get_nbse_mols(dir: str) -> list[Chem.Mol]:
    files = [file for file in listdir(dir) if not file.startswith('benchmark_result') and file.endswith('.sdf')]
    return flatten([sdf_to_smiles(os.path.join(dir, file)) for file in files])


def get_nbse_smiles(dir: str) -> list[str]:
    return [Chem.MolToSmiles(mol) for mol in get_nbse_mols(dir)]


@dataclass
class Result:
    name: str
    took: float
    conformers: int
    local_score: float
    siena_score: float
    avg_conformer_score: float
    def header(self):
        return 'name,took,conformers,local_similarity,avg_conformer_rmsd,siena_rmsd'

    def __str__(self):
        return '{},{},{},{},{},{}'.format(self.name, self.took, self.conformers, self.local_score, self.avg_conformer_score, self.siena_score)


def benchmark_nbse_ensemble(name: str) -> Result:
    directory = os.path.join(nbse_dir, name, 'ligands')
    orig = get_nbse_mols(directory)

    in_file = os.path.join(directory, 'benchmark_input.smi')
    if not os.path.isfile(in_file):
        smiles = get_nbse_smiles(directory)
        with open(os.path.join(directory, 'benchmark_input.smi'), 'w') as f:
            for smi in smiles:
                f.write(smi + '\n')
            f.flush()

    outfile = tempfile.NamedTemporaryFile(mode='r', delete=False, suffix='.sdf')
    cmd_args = [coaler_bin,
                '-f', in_file,
                '-o', outfile.name,
                '-j', str(threads),
                '--conformers', str(num_conformers)]

    start = time.time()

    coaler = subprocess.Popen(cmd_args, stdout=(sys.stdout if verbose else subprocess.DEVNULL))
    coaler.communicate()

    end = time.time()

    out_mol_suppl = Chem.SDMolSupplier(outfile.name)
    out_mols: list[Chem.Mol] = []
    mol: Chem.Mol
    for mol in out_mol_suppl:
        out_mols.append(mol)

    orig_by_smiles: dict[str, Chem.Mol] = {}
    for mol in orig:
        smi = Chem.MolToSmiles(mol).replace('@', '')
        orig_by_smiles[smi] = mol

    out_mols_by_smiles: dict[str, Chem.Mol] = {}
    for mol in out_mols:
        smi = Chem.MolToSmiles(mol).replace('@', '')
        out_mols_by_smiles[smi] = mol

    avg_conformer_score = 0
    for conf in out_mols_by_smiles.values():
        conf_smiles = Chem.MolToSmiles(conf).replace('@', '')

        try:
            orig_conf = orig_by_smiles[conf_smiles]
            avg_conformer_score += AllChem.AlignMol(conf, orig_conf)
        except:
            print("failed to align conformers")
            avg_conformer_score += 1

    avg_conformer_score /= len(orig)

    orig_sorted = sorted(orig_by_smiles.items())
    _, merged_in = orig_sorted[0]
    for _, mol in orig_sorted[1:]:
        merged_in = Chem.CombineMols(merged_in, mol)

    out_sorted = sorted(out_mols_by_smiles.items())
    _, merged_out = out_sorted[0]
    for _, mol in out_sorted[1:]:
        merged_out = Chem.CombineMols(merged_out, mol)

    aligned = AllChem.AlignMol(merged_out, merged_in)

    print("global alignment score: {}".format(aligned))

    result_writer = Chem.SDWriter(os.path.join(directory, 'benchmark_result.sdf'))
    result_writer.write(merged_in)
    result_writer.write(merged_out)
    result_writer.close()

    print('number of molecules in output: {}'.format(len(out_mols)))

    outfile.close()

    return Result(
        name=name,
        took=end - start,
        conformers=num_conformers,
        siena_score=aligned,
        avg_conformer_score=avg_conformer_score,
        local_score=out_mols[0].GetProp("_Score"),
    )

if __name__ == '__main__':
    reinvesitgate = ['2pqk', '1v48', '4ajn', '4hw2', '4ly9']
    done = ['1d0s', '4ajn','1u0z', '2w0v', '4c4f', '3id8', '4c4f', '3id8', '4nb6', '2zsd']
    working = ['1d0s', '4ajn','1u0z', '2w0v', '4c4f', '3id8', '4c4f', '3id8', '4nb6', '2zsd']

    results = []
    for name in working:
        print("running: {}".format(name))
        results.append(benchmark_nbse_ensemble(name))

    results_csv = None
    if not os.path.isfile('benchmark_results.csv'):
        results_csv = open('benchmark_results.csv', 'w')
        results_csv.write(Result.header(results[0]) + '\n')
    else:
        results_csv = open('benchmark_results.csv', 'a')

    for result in results:
        results_csv.write(str(result) + '\n')

    results_csv.close()

    print(results)
