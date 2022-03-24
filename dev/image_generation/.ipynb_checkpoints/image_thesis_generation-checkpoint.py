import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'

plt.style.use('./thesis_style.mplstyle')



color_list = [CB91_Purple, CB91_Green, CB91_Amber, CB91_Pink,
              CB91_Blue, CB91_Violet]
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

def data_image(data, name, title, xlabel, ylabel, legend):
    
    fig, ax = plt.subplots(tight_layout=True, figsize=(8,4))

    for dataplot in data:
        #dataplot[dataplot.columns[0]].plot(ax=ax, kind="line",use_index=True)    
        #dataplot[dataplot.columns[0]].plot(ax=ax, kind="line",use_index=True)    
        
        plt.plot(dataplot.index, dataplot[dataplot.columns[0]])
    
    first = float(data[0].index[0])
    last = float(data[0].index[-1])


    left, right = plt.xlim()
    
    # This is an exception to make the Jitter and Shimmer Frequency figure look nice
    if (last == 24000):
        number_of_ticks = 9
    else:
        number_of_ticks = 8
    
    # The original axis start before the first value and end after the last value. To make it
    # start from 0 to the last value, the range has to go from 0 to right + left
    if (left <= 0):
        xaxis_range = np.linspace(left, right + left, number_of_ticks)
    else:
        xaxis_range = np.linspace(left, right, number_of_ticks)
    xaxis_ticks = np.linspace(0, last, len(xaxis_range))

    # If ticks values are to high, round to the nearest int
    if (xaxis_ticks.max() > 10):
        xaxis_ticks = xaxis_ticks.astype(int)
    else:
        xaxis_ticks = np.around(xaxis_ticks, decimals = 1)
    
    if (left <= 0):
        plt.xlim(left,right + left)
    else:
        plt.xlim(left,right)
    plt.xticks(xaxis_range,xaxis_ticks)
    
    ax.set_title(title, fontsize='large')
    ax.set_xlabel(xlabel, fontsize='medium')
    ax.set_ylabel(ylabel, fontsize='medium')

    ax.legend([legend])#, bbox_to_anchor=(0, -0.4, 1, 0), loc="lower left", mode="expand", ncol=3)

    plt.savefig("thesis/" + name + ".png", bbox_inches='tight', format='png')
    plt.close(fig)
    

original_spectrum = pd.read_csv("../outputs/figures/original_spectrum.csv",index_col=0)
modified_spectrum = pd.read_csv("../outputs/figures/modified_spectrum.csv",index_col=0)
spectrum_modifer = pd.read_csv("../outputs/figures/spectrum_modifier.csv",index_col=0)
rosenberg = pd.read_csv("../outputs/figures/rosenberg_rd.csv",index_col=0)


data_image([original_spectrum], "original_spectrum", "Original Recovered Spectrum", "Frequency [Hz]", "Amplitude [-]","Original Spectrum")
data_image([modified_spectrum], "modified_spectrum", "Modified Recovered Spectrum", "Frequency [Hz]", "Amplitude [-]","Modified Spectrum")
data_image([spectrum_modifer], "spectrum_modifer", "Spectrum Modifier", "Sample [-]", "Amplitude [-]","Modifier")
data_image([rosenberg], "rosenberg_rd", "Rosenberg Pulse RD Parameter Comparison", "Sample [-]", "Amplitude [-]",["Modifier1","Modifier2","Modifie3"])