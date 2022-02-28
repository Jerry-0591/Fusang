# Fusang <img align="right" src="https://github.com/Jerry-0591/Fusang/blob/main/logo.jpg" width="170" height="170"/>
💜 [fusang.cibr.ac.cn](fusang.cibr.ac.cn/) 💜

Fusang is a framework used for the reconstruction of phylogenetic tree via deep learning methods. For current version, it supports the reconstruction of MSA with 4-40 taxas and the length of it should be less than 10,000.

## Hardware requirements 

This repository can be run on both CPU and GPU environment, but we recommend users use GPU for accelerating. More details can be seen from Environment_setting.md

The limit usage of memory is ~24GB for current repository, for most cases, the memory usage is less than 20GB.

## Software requirements

The configuration of the environment see **Environment_setting.md** of this repository

## Example of usage

### 1. Quick start

You can run Fusang using default parameter setting through the command as follows:

```
python fusang.py --msa_dir /path/to/your_msa.fas --save_prefix your_prefix
```

An example of command as follows:

```
python fusang.py --msa_dir ./example_msa/msa1.fas --save_prefix dl_output_1
```

This command will do phylogenetic reconstruction of your MSA file, the result will be saved in file with the prefix that you set in `--save_prefix`

The meaning of these two mandatory parameter:

`--msa_dir` The path to MSA file,  for current version of Fusang, we support both fasta and phylip format of MSA. The example of current MSA format can be seen in the directory of `example_msa`

`--save_prefix`  The prefix of output file, the predicted tree will be saved on the directory of `dl_output` , with the prefix that set in this parameter. You can see `example_dl_output` to find the example of predicted tree.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### 2. Parameter setting

You can set the parameters as follows for specific scenario

`--beam_size` The size of beam in the beam search procedure, default setting is 3. 

`--sequence_type` The size of beam in the beam search procedure  

`--branch_model` The size of beam in the beam search procedure  

`--window coverage` The size of beam in the beam search procedure  

More examples see **Document.md**



## Meaning of each file in this repository

`dl_model` The directory that saves the model of deep learning

`example_dl_output` The directory that saves the example of predicted tree

`example_msa` The directory that saves the example of input msa file

`fusang.py` The code for tree reconstruction 

`requirements.txt` This file will be used for setting environment

`example_command.sh` This command will generate the prdicted tree (in `example_dl_output` ) of input msa file (in `example_msa` ) 

