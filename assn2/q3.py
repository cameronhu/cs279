import math as math
import pandas as pd

kb = 1.38 * 10**-23
dict = {'x': [100, -5], 'y': [10, -10], 'z': [500, -.5]}
temps = [1, 100, 310]

def calc_probability(microstates, pe, temperature):
    exp = (-pe*10**-22) / (kb * temperature)
    return microstates * math.pow(math.e, exp)

def q3():
    prob_dict = {}
    for macrostate in dict.keys():
        macrostate_chars = dict[macrostate]
        prob_dict[macrostate] = []
        for temp in temps:
            this_prob = calc_probability(macrostate_chars[0], macrostate_chars[1], temp)
            print(f"Probability of macrostate {macrostate} at tempearature {temp} is: {this_prob} ")
            prob_dict[macrostate].append(this_prob)
    df = pd.DataFrame(prob_dict, index = temps)
    df.to_csv("probabilities")
    return(df)
