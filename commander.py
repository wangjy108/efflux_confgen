import sys
from pathlib import Path
import subprocess
import time
import os
import pandas as pd

from runner import main as executor

from dp.launching.typing import BaseModel, Field, OutputDirectory, InputFilePath, Int, Float, Enum
from dp.launching.cli import to_runner, default_minimal_exception_handler


class efflux_mediator(str, Enum):
    bcrp = "BCRP"
    pgp = "Pgp"

class Options(BaseModel):
    Input: InputFilePath = Field(description="The input can either be prepared SDF file with 3D coordinates (all-in-one file) \
                                              or SMILES in csv format (all-in-one file); \
                                              if adopt SDF file, please make sure the name of each conform should be unique;\
                                              if use SMILES csv, 1st column should be SMILES (without header),2nd column could be the corresponding NAME (without header), \
                                                  either COMMA, or TAB seperated")
    EffluxMediator: efflux_mediator = Field(render_type="radio", description="efflux mediator, choose 1")
    EnergyWindow: Float = Field(description="energy window applied, unit in kcal/mol, default is 5.0", default=5.0)
    RMSDCutoff: Float = Field(description="rmsd cutoff applied, unit in Angstrom, default is 2.0", default=2.0)
    NposePerConformer: Int = Field(description="saved n_pose for each conformer, default is 1,", default=1)
    
    output_dir: OutputDirectory = Field(description="output dir to save results", default="./outputs")


class GlobalOptions(Options, BaseModel):
    ...

def main(opts: GlobalOptions):

    #print(opts.Input)
    #print(os.getcwd())

    output_dir = Path(opts.output_dir.get_path())
    output_dir.mkdir(parents=True, exist_ok=True)

    gen = executor(input=opts.Input,
                        energy_window=opts.EnergyWindow,
                        rmsd_cutoff=opts.RMSDCutoff,
                        efflux_mediator=opts.EffluxMediator,
                        n_pose_each_conformer=opts.NposePerConformer).go()
    
    
    #prefix = ".".join(opts.Input.split("/")[0].split(".")[:-1])
    mediator = opts.EffluxMediator.upper()

    if not gen:
        return 
    
    with (output_dir / f"FIX4{mediator}.sdf").open("w+") as f:
    #with (output_dir / "output.sdf").open("w+") as f:
        for each in gen[0]:
            f.write(each)
        
    if gen[1]:
        with (output_dir / "ERROR.smi").open("w+") as d:
            for line in gen[1]:
                d.write(line)
    
    if gen[2]:
        with (output_dir / "ISSUE.list").open("w+") as ee:
            for every in gen[2]:
                ee.write(every)

        #df_error.to_csv(output_dir / "ERROR.smi", index=None)


def to_parser():
    return to_runner(
        GlobalOptions,
        main,
        version='0.0.2',
        exception_handler=default_minimal_exception_handler,
    )

if __name__ == '__main__':
    to_parser()(sys.argv[1:])
    

    
