# Environment setting

1. Install [Python 3.7](https://www.python.org/downloads/release/python-370/)

2. Install virtual environment of python

   ```
   pip install virtualenv
   ```

3. Install [Cuda 10.1](https://developer.nvidia.com/cuda-10.1-download-archive-base)(**optional**, only if you want to **use GPU **to accelerate)

   For my testing of 20 taxas with MSA length about 6.6K, the speed of CPU/GPU as follows:
   
   Only CPU (using 12 cores)      ~3.79min 
   Using GPU (NVIDIA TESLA V100)  ~1.01min 
   
   We recommend using GPU for accelerating.
   
4. Use virtual environment of python, it will create a new folder named `fusang_env` in the current directory

   ```
   virtualenv --python=python3.7 fusang_env
   ```

5. Source the virtual environment of python

   ```
   source path/to/fusang_env/bin/activate
   # path/to/fusang_env means that the directory of fusang_env
   ```
   
6. Enter the folder of **this repository**, install python packages

   ```
   cd path/to/version1   # path/to/version1 means that enter the directory of version1
   pip install -r requirements.txt
   
   # This command may take several hours according to the network situation,  
   # we suggest you use command like screen or nohup in linux for installation
   ```

7. Finish the installation of environment

   Note: <big>**every time**</big> before you start the phylogenetic reconstruction, you need to **activate** the virtual environment first, and after you end the phylogenetic reconstruction, you need to **deactivate** it. The command of it as follows:

   ```
   #activate the virtual environment
   source path/to/fusang_env/bin/activate
   
   #deactivate the virtual environment
   deactivate
   ```

   If you have problem about **python virtual environment**, you can see some introductions of it:

   https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/

