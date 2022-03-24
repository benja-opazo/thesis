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

def data_mean(data, measure):
    data_list = []
    for i in data.index:
        data_list.append(data.loc[i].values.mean())
    return pd.DataFrame(data_list, index=data.index, columns=["Mean " + measure])

def smilextract_data(data, measure, datalabel):
    data_list = []

    for i in data.index:
        data_list.append(float(data[datalabel].loc[i]))

    return  pd.DataFrame(data_list, index=data.index, columns=["Mean " + measure])
    
def data_image(data, name, title, xlabel, ylabel):
    
    fig, ax = plt.subplots(tight_layout=True, figsize=(8,4))

    for dataplot in data:
        #dataplot[dataplot.columns[0]].plot(ax=ax, kind="line",use_index=True)    
        #dataplot[dataplot.columns[0]].plot(ax=ax, kind="line",use_index=True)    
        
        plt.plot(dataplot.index, dataplot[dataplot.columns[0]])
    
    first = float(data[0].index[0])
    last = float(data[0].index[-1])


    left, right = plt.xlim()
    
    # This is an exception to make the Jitter and Shimmer Frequency figure look nice
    if (last == 21950):
        last = 22000
        first = 0
        number_of_ticks = 9
    else:
        number_of_ticks = 8
    
    # The original axis start before the first value and end after the last value. To make it
    # start from 0 to the last value, the range has to go from 0 to right + left
    if (left <= 0):
        xaxis_range = np.linspace(0, right + left, number_of_ticks)
    else:
        xaxis_range = np.linspace(left, right, number_of_ticks)
    xaxis_ticks = np.linspace(first, last, len(xaxis_range))

    # If ticks values are to high, round to the nearest int
    if (xaxis_ticks.max() > 10):
        xaxis_ticks = xaxis_ticks.astype(int)
    else:
        xaxis_ticks = np.around(xaxis_ticks, decimals = 1)
    
#    print("Left: " + str(left) + " Right : " + str(right))
#    if (left <= 0):
#        plt.xlim(0,right + left)
#    else:
#        plt.xlim(left,right)
        
    plt.xticks(xaxis_range,xaxis_ticks)
    
    ax.set_title(title, fontsize='large')
    ax.set_xlabel(xlabel, fontsize='medium')
    ax.set_ylabel(ylabel, fontsize='medium')

    if (len(data) == 4):
        ax.legend(['Female Speech', 'Female Sust', 'Male Speech', 'Male Sust'], bbox_to_anchor=(0, -0.4, 1, 0), loc="lower left", mode="expand", ncol=4)
    elif (len(data) == 2):
        ax.legend(['Female Sust', 'Male Sust'], bbox_to_anchor=(0, -0.4, 1, 0), loc="lower left", mode="expand", ncol=4)

    plt.savefig("images/" + name + ".png", bbox_inches='tight', format='png')
    plt.close(fig)
    
    

# Note that independent_variables and measure is a list of a list, where the first element is the name of the files, and the second element is the pretty name
independent_variables = [["f0_multiplier", "F0 Multiplier"],["jitter_amplitude", "Jitter Amplitude"],["jitter_frequency", "Jitter Frequency"],["rd", "Rd Parameter"],["rpp_k", "Pulse Amplitude"],["shimmer_amplitude", "Shimmer Amplitude"],["shimmer_frequency", "Shimmer Frequency"]]
sound_names = ["female_runn", "female_sust", "male_runn", "male_sust"]
objective_measure = [["cpp", "CPP", "dB"], ["ps", "Peak Slope", "-"], ["pesq", "PESQ", "-"], ["se", "Spectral Envelope", "dB"], ["hnr", "Harmonic to Noise Ratio", "-"], ["jitter", "Jitter", "%"], ["shimmer", "Shimmer", "%"]]

file_prefix = "../outputs/cleaned/"

#independent_variables = [independent_variables[3]]
#objective_measure = [objective_measure[0]]


for ind_var in independent_variables:
    for measure in objective_measure:
        
        # This is for the jitter and shimmer exception
        if (len(measure) == 3):
        	measure_aux = measure
        input_data = []
        mean_input_data =[]

        # This for creates a list with input_data of every sound
        for sounds in sound_names:
            print(file_prefix + ind_var[0] + "_" + measure[0] + "_" + sounds + ".csv")
            if (measure[0] == "hnr"):
                input_data.append((pd.read_csv(file_prefix + ind_var[0] + "_" + measure[0] + "_" + sounds + ".csv",index_col=0)))
            elif (measure[0]=="jitter" or measure[0]=="shimmer"):
                input_data.append((pd.read_csv(file_prefix + ind_var[0] + "_" + "jitter_and_shimmer" + "_" + sounds + ".csv",index_col=0)))
            elif (measure[0]=="se"):
                if (sounds == "female_sust" or sounds == "male_sust"):
                    if (ind_var[0] == "rd" or ind_var[0] == "rpp_k" or ind_var[0] == "f0_multiplier"):
                        input_data.append((pd.read_csv(file_prefix + ind_var[0] + "_" + measure[0] + "_" + sounds + ".csv",index_col=0).transpose()))    
            else:
                input_data.append((pd.read_csv(file_prefix + ind_var[0] + "_" + measure[0] + "_" + sounds + ".csv",index_col=0).transpose()))

        # This for creates a list with the mean of each data, ready to generate the figure
        for data in input_data:
            if (measure[0] == "hnr"):
                mean_input_data.append(smilextract_data(data, measure[1], 'HNR_sma_amean'))
            elif (measure[0]=="jitter"):
                mean_input_data.append(smilextract_data(data, measure[1], 'jitterLocal_sma_amean'))
            elif (measure[0]=="shimmer"):
                mean_input_data.append(smilextract_data(data, measure[1], 'shimmerLocal_sma_amean'))
            else:
                mean_input_data.append(data_mean(data, measure[1]))

        # Figure generation
        if len(mean_input_data) > 0:
            data_image(mean_input_data , ind_var[0] + "_" + measure[0] , ind_var[1] +  " vs " + measure[1], ind_var[1] , "Mean " + measure[1] + " [" + measure[2] + "]")
