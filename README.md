# Biological Sequence Learning with Conserved Domain
This project contains deep learning experiments on biological sequences (amino acid or nucleotide sequences) using conserved domain information. 
Please refer to the [DESCRIPTION](docs/DESCRIPTION.md) file for more technical details.


## Table of Contents  
-   [Getting Started](#getting-started)
-   [Template Layout](#template-layout)
-   [Future Tasks](#future-tasks)
-   [Authors](#authors)
-   [License](#license)


---
## Getting Started
This project uses [Poetry](https://python-poetry.org/) for dependency management. 
Check 'pyproject.toml' for auto-generated dependencies. 
```bash
conda install poetry  # or 'pip install poetry'
poetry install
```


---
## Template Layout
The project layout with the usage for each folder is shown below:
```text
dl-project-template
.
|
├── LICENSE.md
├── README.md
├── DESCRIPTION.md      # project and experiment description
├── makefile            # makefile for various commands (install, train, pytest, mypy, lint, etc.) 
├── mypy.ini            # MyPy type checking configurations
├── pylint.rc           # Pylint code quality checking configurations
├── pyproject.toml      # Poetry project and environment configurations
|
├── data
|   ├── ...             # data reference files (index, readme, etc.)
│   ├── raw             # untreated data directly downloaded from source
│   ├── interim         # intermediate data processing results
│   └── processed       # processed data (features and targets) ready for learning
|
├── notebooks           # Jupyter Notebooks for data processing experiments and visualization
├── scripts             # shell scripts for data download and processing
│── src    
│   │── ...             # top-level scripts for training, testing and global tasks
│   ├── configs         # configuration files for deep learning experiments
│   ├── datasets        # dataset classes and other processing functions
│   ├── modules         # activations, layers, modules, and networks (subclass of torch.nn.Module)
│   ├── optimization    # deep learning optimizers and schedulers
│   └── utilities       # other useful functions and classes
├── tests               # unit tests module for ./src
│
├── docs                # documentation files (*.txt, *.doc, *.jpeg, etc.)
├── logs                # logs for deep learning experiments
└── models              # saved models with optimizer states
```


---
## Future Tasks
<!-- TODO -->


---
## Authors
* Xiaotian Duan (Email: xduan7 at uchicago.edu)


---
## License
This project is licensed under the MIT License - check the [LICENSE](LICENSE.md) file for more details.

