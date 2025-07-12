Code repository for:  
F. M. Barbero, T. Lenc, N. Jacoby, R. Polak, M. Varlet, S. Nozaradan, Revealing rhythm categorization in human brain activity. *Science Advances*, doi: 10.1126/sciadv.adu9838 (2025).

After downloading the data, go to `get_par.m` and update the variable `data_path` so it points to the base folder where your data is. 
## Dependencies
- [letswave6](https://github.com/NOCIONS/letswave6)
- [letswave7](https://github.com/NOCIONS/letswave7)
- [rnb_tools (commit `cef0119`)](https://github.com/TomasLenc/rnb_tools/tree/cef0119ac0fbfa442e0e41ef8c6b10e1addac5c3)

After downloading each library, make sure to update the path to it in `get_par.m`.  

## Pipeline

Make sure your matlab working directory is set to the root folder of the repo (where `get_par.m` file is located). 


### Auditory nerve model 

First you need to compile some mex files. Go to `lib/urear` and follow the README.  

When done, execute scripts in the `auditory_nerve_model` subfolder. 

### Preprocessing

Execute scripts in the `preprocessing` subfolder.  

Each script is run separately for each individual subject (there are some manual steps involved).  


### Main analyses

After preprocessing all subjects, run the scripts in the subfolder `summary_scripts`. 

