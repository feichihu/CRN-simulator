import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('res.csv')

def plots(data):
    fig = plt.figure()
    plt.title("Simulation - CSV")
    plt.xlabel("time (seconds)")
    plt.ylabel("molecule counts")
    axi = []
    names = data.columns[2:]
    l = []
    for i in range(len(names)-1):
        axi.append(fig.add_axes())
    for i in range(len(names)-1):
        axi[i] = plt.plot(data['time'], data[names[i+1]])
        l.append((names[i+1]))
    plt.legend(l)
    print(names)
    return True
    
plots(data)
plt.show()
