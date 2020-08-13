# Deep Learning pn Biological Sequences
This project contains deep learning experiments, mostly unsupervised learning, on biological sequences (protein sequences and gene sequences). 


## Table of Contents  
- [Dependencies](#dependencies)
- [Getting Started](#getting-started)
- [Template Layout](#template-layout)
- [Future Tasks](#future-tasks)
- [Authors](#authors)
- [License](#license)


## Dependencies
```text
required:
    python>=3.8
    numpy>=1.18
    pandas>=1.0
    torch>=1.4
    scikit-learn>=0.22

optional:
    poetry>=0.12
    flake8>=3.7
    pylint>=2.4
    mypy>=0.76
    pytest>=5.3
    GPUtil>=1.4
```


## Getting Started


## Template Layout
The project layout with the usage for each folder is shown below:
```text
dl-project-template
.
|
├── LICENSE.md
├── README.md
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
├── notebooks           # Jupyter Notebooks (mostly for data processing and visualization)
│── src    
│   │── ...             # top-level scripts for training, testing and downloading
│   ├── configs         # configuration files for deep learning experiments
│   ├── data_processes  # data processing functions and classes (cleaning, validation, imputation etc.)
│   ├── modules         # activations, layers, modules, and networks (subclass of torch.nn.Module)
│   ├── optimization    # deep learning optimizers and schedulers
│   └── utilities       # other useful functions and classes
├── tests               # unit tests module for ./src
│
├── docs                # documentation files (*.txt, *.doc, *.jpeg, etc.)
├── logs                # logs for deep learning experiments
└── models              # saved models with optimizer states
```

## Future Tasks
- [ ] ML/DL projects process flowchart 
    - [ ] definition of several major steps
    - [ ] clarify motivation and deliverables
- [ ] small example for demonstration (omniglot?)


## Authors
* Xiaotian Duan (Email: xduan7 at uchicago.edu)


## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for more details.

