import os
from copy import deepcopy
import numpy as np
from opfunu.cec_based import cec2020
from pyDOE2 import lhs


PopSize = 100
DimSize = 10
LB = [-100] * DimSize
UB = [100] * DimSize
TrialRuns = 30
MaxFEs = DimSize * 10000
curFEs = 0

Epsilon = 2 ** (-52)
MaxIter = int(MaxFEs / PopSize)
curIter = 0

Pop = np.zeros((PopSize, DimSize))
Off = np.zeros((PopSize, DimSize))

FitPop = np.zeros(PopSize)
FitOff = np.zeros(PopSize)

Func_num = 0
SuiteName = "CEC2020"

X_best = None
X_best_Fit = float("inf")


def Distance(indiX, indiY):
    global DimSize
    dis = 0
    for i in range(DimSize):
        dis += abs(indiX[i] - indiY[i])
    return dis + Epsilon


# initialize the Pop randomly
def Initialization(func):
    global Pop, FitPop, curFEs, DimSize, X_best, X_best_Fit
    Pop = lhs(DimSize, samples=PopSize)
    for i in range(PopSize):
        for j in range(DimSize):
            Pop[i][j] = LB[j] + (UB[j] - LB[j]) * Pop[i][j]
        FitPop[i] = func.evaluate(Pop[i])
        curFEs += 1
    X_best_Fit = min(FitPop)
    X_best = deepcopy(Pop[np.argmin(FitPop)])



def RIME(func):
    global Pop, FitPop, Off, FitOff, curIter, MaxIter, LB, UB, PopSize, DimSize, curFEs, X_best, X_best_Fit
    w = 5
    Off = deepcopy(Pop)
    rime_factor = np.random.uniform(-1, 1) * np.cos(np.pi * (curIter + 1) / (MaxIter / 10)) * (
            1 - np.round((curIter + 1) * w / MaxIter) / w)
    E = np.sqrt((curIter + 1) / MaxIter)
    NorFit = FitPop / np.linalg.norm(FitPop, axis=0, keepdims=True)

    for i in range(PopSize):
        r1, r2 = np.random.choice(list(range(PopSize)), 2, replace=False)
        for j in range(DimSize):
            if np.random.rand() < E:  # Soft rime
                Off[i][j] = X_best[j] + rime_factor * (np.random.rand() * (UB[j] - LB[j]) + LB[j])
            if np.random.rand() < NorFit[i]:  # Hard rime
                Off[i][j] = X_best[j] + NorFit[i] * (Pop[r1][j] - Pop[r2][j])
        Off[i] = np.clip(Off[i], LB, UB)
        FitOff[i] = func.evaluate(Off[i])
        curFEs += 1
        if FitOff[i] < FitPop[i]:
            FitPop[i] = deepcopy(FitOff[i])
            Pop[i] = deepcopy(Off[i])
            if FitOff[i] < X_best_Fit:
                X_best_Fit = deepcopy(FitOff[i])
                X_best = deepcopy(Off[i])
        else:
            dis = Distance(Pop[i], Off[i])
            delta = abs(FitOff[i] - FitPop[i])
            threshold = np.exp(-(delta / dis))
            if np.random.rand() < threshold:
                FitPop[i] = deepcopy(FitOff[i])
                Pop[i] = deepcopy(Off[i])


def RunRIME(func):
    global curFEs, curIter, MaxFEs, TrialRuns, Pop, FitPop, DimSize
    All_Trial_Best = []
    for i in range(TrialRuns):
        Best_list = []
        curFEs = 0
        curIter = 0
        Initialization(func)
        Best_list.append(X_best_Fit)
        np.random.seed(2022 + 88 * i)
        while curFEs < MaxFEs:
            RIME(func)
            curIter += 1
            Best_list.append(X_best_Fit)
        All_Trial_Best.append(Best_list)
    np.savetxt("./" + SuiteName + "_Data/F" + str(Func_num) + "_" + str(DimSize) + "D.csv", All_Trial_Best,
               delimiter=",")


def main(dim):
    global Func_num, DimSize, Pop, MaxFEs, SuiteName, LB, UB
    DimSize = dim
    Pop = np.zeros((PopSize, dim))
    MaxFEs = dim * 10000
    LB = [-100] * dim
    UB = [100] * dim

    CEC2020Funcs = [cec2020.F12020(dim), cec2020.F22020(dim), cec2020.F32020(dim), cec2020.F42020(dim),
                    cec2020.F52020(dim), cec2020.F62020(dim), cec2020.F72020(dim), cec2020.F82020(dim),
                    cec2020.F92020(dim), cec2020.F102020(dim)]
    CEC_dict = {"CEC2020": CEC2020Funcs}
    for suite in CEC_dict.keys():
        Func_num = 0
        SuiteName = suite
        for i in range(len(CEC2020Funcs)):
            Func_num = i + 1
            RunRIME(CEC_dict[SuiteName][i])


if __name__ == "__main__":
    if os.path.exists('./CEC2020_Data') == False:
        os.makedirs('./CEC2020_Data')
    Dims = [10, 30, 50]
    for Dim in Dims:
        main(Dim)
