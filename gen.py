#!/usr/bin/env python
# coding: utf-8

"""
// author: Wang (Max) Jiayue 
// email: wangjy108@outlook.com
"""

"""
current version follow the strain workflow strategy
"""

import argparse
import logging
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import os
import pandas as pd
import numpy as np
import subprocess
from multiprocessing import cpu_count
from util.ConfGenbyMM import ConfGen
#from util.Cluster import cluster
from util.OptbySQM import System as sysopt
#from util.SPcalc import System as syssp
from util.Align import Align as align
from util.ConfRelaxbySQM import System as MDsample
from util.CalRMSD import RMSD as RMSD
import shutil

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)


class main():
    def __init__(self, **args):
        #self.main_dir = os.getcwd()
        self.input_smi = args["input_smi"]  # str   
        self.input_mol_label = args["input_mol_label"]

        self.n_thread = cpu_count()

        try:
            self.energy_window = float(args["energy_window"])
        except Exception as e:
            self.energy_window = 5.0
        
        try:
            self.rmsd_cutoff = float(args["rmsd_cutoff"])
        except Exception as e:
            self.rmsd_cutoff = 2.0
    
    def gen_initial_conformer_by_rdkitMM(self):
        ConfGen(input_smi_file=self.input_smi, \
                method="MMFF94", \
                genConfNum=20, \
                saveConfNum=1, \
                fileName=f"_INITIAL_{self.input_mol_label}.sdf").run()
        
        if os.path.isfile(f"_INITIAL_{self.input_mol_label}.sdf") \
            and os.path.getsize(f"_INITIAL_{self.input_mol_label}.sdf"):
            return f"_INITIAL_{self.input_mol_label}.sdf"
        else:
            return None             
    
    def gen_initial_conformer_by_obabel(self):
        gen_cmd = f'obabel -:{self.input_smi} -O _INITIAL_{self.input_mol_label}.sdf --gen3D'
        #(_, _) = subprocess.getstatusoutput(gen_cmd)
        try:
            p = subprocess.run(gen_cmd.split(), timeout=20, check=True, capture_output=True)
        except subprocess.TimeoutExpired:
            logging.info("Timeout with obabel gen")
        
        if os.path.isfile(f"_INITIAL_{self.input_mol_label}.sdf") \
            and os.path.getsize(f"_INITIAL_{self.input_mol_label}.sdf"):
            return f"_INITIAL_{self.input_mol_label}.sdf"
        else:
            return None

    def intial_conformer(self):
        try:
            get_initial = self.gen_initial_conformer_by_rdkitMM()
        except Exception as e:
            get_initial = None
        if not get_initial:
            get_initial = self.gen_initial_conformer_by_obabel()

        return get_initial
    
    def sample(self, input_sdf):

        ## initial
        _dic_rotableBond = {0: 100, 1: 200}
        try:
            get_mol = [mm for mm in Chem.SDMolSupplier(input_sdf) if mm][0]
        except Exception as e:
            get_mol = None
            rotable_bond_index = 0
            logging.info("Failed at generate initial conformer, abort")
            return None
        else:
            rotable_bond = rdMolDescriptors.CalcNumRotatableBonds(get_mol)
            _def_func = lambda x: 0 if max(5, x) == 5 else 1
            rotable_bond_index = _def_func(rotable_bond)
        
        N_conformer = _dic_rotableBond[rotable_bond_index]

        ## perform initial optimization to relax MM bad contact
        try:
            get_initial_opt = sysopt(input_sdf=input_sdf, 
                                    HA_constrain=True,
                                    if_write_sdf=True).run()
        except Exception as e:
            logging.info("Failed at initial optimization after conformer generation, abort")
            return None
        
        if os.path.isfile("_OPT.sdf") and os.path.getsize("_OPT.sdf"):
            optimized_input_sdf = f"{self.input_mol_label}.initial_opt.sdf"
            os.system(f"mv _OPT.sdf {optimized_input_sdf}")
            get_charge = get_initial_opt[0].GetProp("charge")
        else:
            logging.info("Failed at initial optimization after conformer generation, abort")
            return None
        
        ## file name saved with "SAVE", should be divided in different folder if multirun

        MDsample(input_sdf=optimized_input_sdf, 
                 save_frame=N_conformer, 
                 define_charge = get_charge).run()
        
        if not (os.path.isfile("SAVE.sdf") and os.path.getsize("SAVE.sdf")):
            logging.info("Failed at MD sampling after initial optimization, abort")
            return None
        
        opted = sysopt(input_sdf="SAVE.sdf", 
                       rmsd_cutoff=self.rmsd_cutoff,
                        # save_n=self.sp_max_n,
                       HA_constrain=True,
                       define_charge=get_charge,
                       if_write_sdf=False).run()
        
        if not opted:
            logging.info("Failed at post-sampling optimization")
            return None

        sorted_opt = sorted(opted, key=lambda x: float(x.GetProp("Energy_xtb")))
        
        i = 1

        saved = sorted_opt[:1]

        while i < len(sorted_opt):
            dif_ee = abs(float(sorted_opt[i].GetProp("Energy_xtb")) - float(saved[-1].GetProp("Energy_xtb"))) * 627.51
            if dif_ee > self.energy_window:
                saved.append(sorted_opt[i])

            i += 1
        #saved = sorted_opt[:1]

        cc = Chem.SDWriter(f"{self.input_mol_label}.sdf")
        for idx, each in enumerate(saved):
            each.SetProp("_Name", f"{self.input_mol_label}_{idx}")
            cc.write(each)
        cc.close()

        #logging.info(f"Generated mol saved in {self.input_mol_label}.sdf")
        return f"{self.input_mol_label}.sdf"
        

        """
        _dic = {f"{self.input_mol_label}_0": optimized_poses[0]}
        
        # align and get_rmsd to reduce redundant
        rmsd_matrix = np.full((len(optimized_poses), len(optimized_poses)), self.rmsd_cutoff)

        i = 0
        base = optimized_poses[0]

        while i < len(optimized_poses):
            target = optimized_poses[i]
            j = i + 1
            while j < len(optimized_poses):
                search = optimized_poses[j]
                get_search_aligned = align(SearchMolObj=search,
                                             RefMolObj=target,
                                             method="crippen3D").run()
                get_rmsd = RMSD(rdmolobj_mol=get_search_aligned,
                                rdmolobj_ref=target,
                                method="selfWhole").run()
                rmsd_matrix[i, j] = get_rmsd

                j += 1    
            i += 1
        
        x, y = np.where(rmsd_matrix < self.rmsd_cutoff)

        for i in range(1, len(optimized_poses)):
            if i not in list(y):
                aligned_this = align(SearchMolObj=optimized_poses[i],
                                             RefMolObj=base,
                                             method="crippen3D").run()
                _dic.setdefault(f"{self.input_mol_label}_{i}", aligned_this)
        
        cc = Chem.SDWriter(f"{self.input_mol_label}.sdf")
        for kk, vv in _dic.items():
            vv.SetProp("_Name", kk)
            cc.write(vv)
        cc.close()

        logging.info(f"Generated mol saved in {self.input_mol_label}.sdf")
        return f"{self.input_mol_label}.sdf"
        """
    
    def run(self):
        ## initial_gen
        initial_conformer = self.intial_conformer()
        gen_confomer_filename = self.sample(initial_conformer)

        return gen_confomer_filename

def gen_former(input_file,
               energy_window,
               rmsd_cutoff):
    main_dir = os.getcwd()
    df = pd.read_csv(input_file, sep="\\s+", header=None)
    df_header = pd.read_csv(input_file, sep="\\s+")

    input_file_name = input_file.split(".")[0]
    
    if list(df_header.columns)[0] == df.iloc[:,0].to_list()[0]:
        in_use_df = df
    else:
        in_use_df = df_header
    
    smile_list = in_use_df.iloc[:, 0]
    try:
        label_list = in_use_df.iloc[:, 1]
    except Exception as e:
        label_list = [f"{input_file_name}_{ii}" for ii in range(len(smile_list))]

    _dic_error = {}

    #root_save = os.path.join(main_dir, f"{input_file_name}.sdf")

    #if not os.path.isfile(root_save):
    #    os.system(f"touch {root_save}")

    gen_sdf_content = []

    for idx, smi in enumerate(smile_list):
        work_dir = os.path.join(main_dir, label_list[idx])
        if not os.path.exists(work_dir):
            os.mkdir(work_dir) 

        os.chdir(work_dir)
        logging.info(f"working with {label_list[idx]}")

        saved_sdf_filename = main(input_smi=smi,
             input_mol_label=label_list[idx],
             energy_window=energy_window,
             rmsd_cutoff=rmsd_cutoff).run() 

        if saved_sdf_filename:
            with open(saved_sdf_filename, "r+") as cc:
                gen_sdf_content += [ll for ll in cc.readlines()]
            
            #cmd = f"cat {saved_sdf_filename} >> {root_save}"
            #(_, _) = subprocess.getstatusoutput(cmd)
            
        else:
            _dic_error.setdefault(label_list[idx], smi)

        os.chdir(main_dir)
        shutil.rmtree(f"{work_dir}")
    
    #if _dic_error:
        #df = pd.DataFrame({"smiles": [vv for vv in _dic_error.values()],
        #                   "label": [kk for kk in _dic_error.keys()]})
        #df.to_csv("ERROR.smi", index=None, sep="\t")
    
    #logging.info("Finish, check *.sdf for generated poses, check ERROR.smi for failed smiles if exists")
    return gen_sdf_content, _dic_error
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run conf gen')
    parser.add_argument('--input_file', type=str, required=True)
    parser.add_argument('--energy_window', type=float, default=5.0)
    parser.add_argument('--rmsd_cutoff', type=float, default=2.0)

    args = parser.parse_args()

    gen_sdf_content, _dic_error = gen_former(input_file=args.input_file,
                                                energy_window=args.energy_window,
                                                rmsd_cutoff=args.rmsd_cutoff)

    with open("Output.sdf", "w+") as c:
        for each in gen_sdf_content:
            c.write(each)
    
    if _dic_error:
        with open("Error.smi", "w+") as ee:
            for kk, vv in _dic_error.items():
                ee.write(f"{kk}\t{vv}\n")

            
