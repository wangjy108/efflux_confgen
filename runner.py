import os
import logging
from rdkit import Chem
import numpy as np
import pandas as pd
from gen import gen_former
from dock import main as sdf2pdb
from pdb2sdf import main as pdb2sdf

import argparse
import configparser

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self, **args):
        self.input = args["input"]
        self.energy_window = args["energy_window"]
        self.rmsd_cutoff = args["rmsd_cutoff"]

        self.receptor = args["efflux_mediator"].lower()
        self.max_pose_n = args["n_pose_each_conformer"]

        self.prefix = ".".join(self.input.split("/")[-1].split(".")[:-1])
    
    
    def go(self):
        main_dir = os.getcwd()

        real_input_smi = ""
        sdf_4_docking = []
        ## try if 3D coordinates involved 
        if self.input.split(".")[-1].lower() == "sdf":
            input_as_smi = {}
            try:
                mol_set = [mm for mm in Chem.SDMolSupplier(self.input, removeHs=False) if mm]
            except Exception as e:
                logging.info("Bad SDF input, abort")
                os.chdir(main_dir)
                return None
            
            for idx, eachMol in enumerate(mol_set):
                try:
                    get_name = eachMol.GetProp("_Name")
                except Exception as e:
                    get_name = f"{self.prefix}_{idx}"
                if not get_name:
                    get_name = f"{self.prefix}_{idx}"
                
                try:
                    input_as_smi[get_name]
                except Exception as e:
                    pass
                else:
                    get_name = f"{self.prefix}_{idx}"
                
                mol_xyz_z = eachMol.GetConformer().GetPositions()[:, 2]

                if abs(np.std(mol_xyz_z) - 0.0) < 1e-3:
                    input_as_smi.setdefault(get_name, Chem.MolToSmiles(eachMol))
                else:
                    sdf_4_docking.append(eachMol)
            
            if input_as_smi:
                real_input_smi = f"TEMP_{self.prefix}.smi"

                df = pd.DataFrame({"SMILES": [vv for vv in input_as_smi.values()],
                                   "NAME": [kk for kk in input_as_smi.keys()]})
                
                df.to_csv(real_input_smi, header=None, index=None, sep=" ")

        else:
            real_input_smi = self.input
        
        error_from_gensmi = []
        

        if real_input_smi:
            saved_sdf_filename = None
            logging.info("----> Step 1: Transforming 1D SMILES to 3D coordinates")

            gen_sdf_content, _dic_error = gen_former(input_file = real_input_smi,
                                            energy_window = self.energy_window,
                                            rmsd_cutoff = self.rmsd_cutoff)

            if gen_sdf_content:
                saved_sdf_filename = f"TEMP_{self.prefix}.sdf"
                with open(saved_sdf_filename, "w+") as cc:
                    for line in gen_sdf_content:
                        cc.write(line)
                if _dic_error:
                    for kk, vv in _dic_error.items():
                        error_from_gensmi.append(f"{kk}\t{vv}\n")

                logging.info("----> Step 1: Done")
            else:
                logging.info("----> Step 1: Error")
            #    return None
            
            if sdf_4_docking:
                if saved_sdf_filename:
                    ready = open(saved_sdf_filename, "a+")
                    cc = Chem.SDWriter(ready)
                    for each in sdf_4_docking:
                        cc.write(each)
                    
                    cc.close()
                    ready.close()
                    #os.system(f"mv {saved_sdf_filename} TEMP_{self.prefix}.sdf")

                else:
                    cc = Chem.SDWriter(f"TEMP_{self.prefix}.sdf")
                    for each in sdf_4_docking:
                        cc.write(each)
                    
                    cc.close()

        else:
            if not sdf_4_docking:
                logging.info("Nothing to do, abort")
                os.chdir(main_dir)
                return None
            
            cc = Chem.SDWriter(f"TEMP_{self.prefix}.sdf")
            for each in sdf_4_docking:
                cc.write(each)
            
            cc.close()

            #saved_sdf_filename = f"_TEMP_{self.prefix}.sdf"
            
            logging.info("Skipping Step 1 as input mols have 3D coordinates")
        
        logging.info(f"----> Step 2: Fixing 3D coordicates into unipfied axis for {self.receptor.upper()}")

        mediator = {
            "bcrp": "/root/commander/config/6vxj_BCRP.pdb",
            "pgp": "/root/commander/config/6qex_Pgp.pdb"
        }

        try:
            receptor_name = mediator[self.receptor]
        except Exception as e:
            logging.info("----> Wrong mediator, choose from [bcrp, Pgp]")
            logging.info("----> Step 2: Error, abort")
            os.chdir(main_dir)
            return None
    

        sdf2pdb(receptor_name=receptor_name,
                ligand_db_name=f"TEMP_{self.prefix}.sdf",
                rmsd_cutoff=self.rmsd_cutoff,
                max_n_poses=self.max_pose_n).run()
        
        if os.path.exists("DockingPose"):
            logging.info("----> Step 2: Done")

        else:
            logging.info("----> Fix failed")
            logging.info("----> Step 2: Error, abort")
            os.chdir(main_dir)
            return None
        
        ## shifting and prepare
        finale_sdf, issued = pdb2sdf(ligand_db_name=f"TEMP_{self.prefix}.sdf").run()

        ## change name_tag
        #name_tag = f"FIX4{self.receptor.upper()}{finale_sdf}"

        #with open(finale_sdf, "r+") as cc:
        #    gen_content = [line for line in cc.readlines()]
            
        os.system("rm -f *TEMP_*")

        if error_from_gensmi:
            logging.info(f"----> Warning: some of the molecules failed at initial preparation, check ERROR.smi")
        
        df_error = []

        if issued:
            logging.info(f"----> Warning: some of the molecules failed at docking, check ISSUE.list")
            for each in issued:
                df_error += f"{each}\n"
                
            #df_error = pd.read_csv("ERROR.smi", header=0)
            #with open("ERROR.smi", "r+") as ff:
            #    df_error = [ll for ll in ff.readlines()]
                
            #os.system("rm -f ERROR.smi")

        if finale_sdf:
            logging.info(f"----> Finish")
            os.chdir(main_dir)
                
            return (finale_sdf, error_from_gensmi, df_error)

        else:
            logging.info(f"----> Something Wrong, check process logging info, and run again")
            os.chdir(main_dir)
            return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='save your ass')
    parser.add_argument('--input', type=str, required=True)
    parser.add_argument('--energy_window', type=float, required=True)
    parser.add_argument('--rmsd_cutoff', type=float, required=True)
    parser.add_argument('--efflux_mediator', type=str, required=True)
    parser.add_argument('--n_pose_each_conformer', type=str, required=True)


    args = parser.parse_args()
    
    cc = main(input=args.input,
         energy_window=args.energy_window,
         rmsd_cutoff=args.rmsd_cutoff,
         efflux_mediator=args.efflux_mediator,
         n_pose_each_conformer=args.n_pose_each_conformer).go()

    if cc:
        with open("output.sdf", "w+") as ww:
            for each in cc[0]:
                ww.write(each)
        
        if cc[1]:
            with open("ERROR.smi", "w+") as e:
                for line in cc[1]:
                    e.write(line)
        
        if cc[2]:
            with open("ISSUE.list", "w+") as ee:
                for every in cc[2]:
                    ee.write(every)
        
            



            
        


        

        



        
        