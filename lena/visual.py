from matplotlib import pyplot as plt

def draw_func(x, f, name):
    
    fig, ax = plt.subplots()

    ax.plot(x, f)
    
    ax.set_xlabel('t')
    ax.set_ylabel(name)
    
    plt.show()
