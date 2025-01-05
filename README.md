# SRIME
SRIME: A strengthened RIME with Latin hypercube sampling and embedded distance-based selection for engineer optimization problems

## Abstract
This paper proposes a strengthened RIME algorithm to tackle continuous optimization problems. RIME is a newly proposed physical-based evolutionary algorithm (EA) inspired by the soft and hard rime growth process of rime-ice, which has a powerful exploitation ability. But in complex optimization problems, RIME will easily trap into local optima and the optimization will become stagnation. Noticing this issue, we introduce three techniques to the original RIME: (1) Latin hypercube sampling replaces the random generator as the initialization strategy, (2) modified hard rime search strategy, and (3) embedded distance-based selection mechanism. We evaluate our proposed SRIME in 10-D, 30-D, 50-D, and 100-D CEC2020 benchmark functions and eight real-world engineering optimization problems with nine state-of-the-art EAs. Experimental and statistical results show that the introduction of three techniques can significantly accelerate the optimization of the RIME algorithm, and SRIME is a competitive optimization technique in real-world applications. Ablation experiments are also provided to analyze our proposed three techniques independently, and the embedded distance-based selection contributes most to the improvement of SRIME. The source code of SRIME can be found in https://github.com/RuiZhong961230/SRIME.

## Citation
@article{Zhong:24,  
author = {Rui, Zhong and Jun, Yu and Chao, Zhang and Masaharu, Munetomo},  
year = {2024},  
month = {01},  
pages = {6721–6740},  
title = {SRIME: A strengthened RIME with Latin hypercube sampling and embedded distance-based selection for engineering optimization problems},  
volume = {36},  
journal = {Neural Computing and Applications},  
doi = {https://doi.org/10.1007/s00521-024-09424-4 },  
}

## Datasets and Libraries
CEC benchmarks are provided by the opfunu library and engineering problems are provided by the enoppy library.


