
def cpp_cleaner(base_name, new_name, index_name, param_min, param_leap, param_max, param_multiplier):
	base_file = open(base_name, 'r')
	new_file = open(new_name, 'w')

	# First we add the indexes
	new_file.write(index_name)

	# param_min, max and leap are used to define the range of the independent variable
	# param_multiplier is used so that range is composed of only integers number
	# The if is because matlab's implementation of range is different
	if (param_max-param_min)/param_leap %1 == 0:
		for i in range(int(param_min*param_multiplier),int((param_max + param_leap )*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))
	else:
		for i in range(int(param_min*param_multiplier),int((param_max)*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))

	new_file.write("\n")

	# Then we add the measurement name
	new_file.write("CPP,")

	# Then we add the values
	lines = base_file.readlines()

	for line in lines:
		new_file.write(line)

	new_file.close()
	base_file.close()


def pesq_cleaner(base_name, new_name, index_name, param_min, param_leap, param_max, param_multiplier):
	base_file = open(base_name, 'r')
	new_file = open(new_name, 'w')

	# First we add the indexes
	new_file.write(index_name)

	# param_min, max and leap are used to define the range of the independent variable
	# param_multiplier is used so that range is composed of only integers number
	# The if is because matlab's implementation of range is different
	if (param_max-param_min)/param_leap %1 == 0:
		for i in range(int(param_min*param_multiplier),int((param_max + param_leap )*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))
	else:
		for i in range(int(param_min*param_multiplier),int((param_max)*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))

	new_file.write("\n")

	# Then we add the measurement name
	new_file.write("PESQ,")

	# Then we add the values
	lines = base_file.readlines()
	for line in lines:
		new_file.write(line)

	new_file.close()
	base_file.close()

def se_cleaner(base_name, new_name, index_name, param_min, param_leap, param_max, param_multiplier):
	base_file = open(base_name, 'r')
	new_file = open(new_name, 'w')

	# First we add the indexes
	new_file.write(index_name)

	# param_min, max and leap are used to define the range of the independent variable
	# param_multiplier is used so that range is composed of only integers number
	# The if is because matlab's implementation of range is different
	if (param_max-param_min)/param_leap %1 == 0:
		for i in range(int(param_min*param_multiplier),int((param_max + param_leap )*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))
	else:
		for i in range(int(param_min*param_multiplier),int((param_max)*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))

	new_file.write("\n")

	# Then we add the measurement name
	new_file.write("SE,")

	# Then we add the values
	lines = base_file.readlines()
	for line in lines:
		new_file.write(line)

	new_file.close()
	base_file.close()

def ps_cleaner(base_name, new_name, index_name, param_min, param_leap, param_max, param_multiplier):
	base_file = open(base_name, 'r')
	new_file = open(new_name, 'w')

	# First we add the indexes
	new_file.write(index_name)


	# param_min, max and leap are used to define the range of the independent variable
	# param_multiplier is used so that range is composed of only integers number
	# The if is because matlab's implementation of range is different
	if (param_max-param_min)/param_leap %1 == 0:
		for i in range(int(param_min*param_multiplier),int((param_max + param_leap )*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))
	else:
		for i in range(int(param_min*param_multiplier),int((param_max)*param_multiplier),int(param_leap*param_multiplier)):
			new_file.write("," + str(i/param_multiplier))

	new_file.write("\n")

	# Then we add the measurement name
	new_file.write("PS,")

	# Then we add the values
	lines = base_file.readlines()

	for line in lines:
		new_file.write(line)

	new_file.close()
	base_file.close()


def hnr_jys_cleaner(base_name, new_name, index_name, param_min, param_leap, param_max, param_multiplier, data_leaps):
	base_file = open(base_name, 'r')
	new_file = open(new_name, 'w')

	# In this case we do the 2 things at the same time

	lines = base_file.readlines()
	line_iter = iter(lines)

	# We add the index
	new_file.write(index_name + ',')
	# We write the line
	new_file.write(next(line_iter))

	#Blank lines
	next(line_iter)

	# The if is because matlab's implementation of range is different
	if (param_max-param_min)/param_leap %1 == 0:
		param_max = param_max + param_leap

	for i in range(int(param_min*param_multiplier),int((param_max)*param_multiplier),int(param_leap*param_multiplier)):

		line_split = []
		line_splitted = []

		individual_value = 0
		average_value_list = []
		#print(i)
		for j in range(data_leaps):
			#print(j)

			line = next(line_iter)
			line_split = line.split(",")
			line_splitted.append(line_split[1:-1]) # doesnt include the 'unknown' and 'class' columns

		for k in range(len(line_splitted[0])):
			individual_value = 0
			for l in range(data_leaps):
				individual_value = individual_value + float(line_splitted[l][k])
			#print(individual_value)
			average_value_list.append(individual_value/data_leaps)

		new_file.write(str(i/param_multiplier) + ',')
		new_file.write('unknown,')
		for writing_line in average_value_list:
			new_file.write(str(writing_line) + ",")
		new_file.write('?\n')


	new_file.close()
	base_file.close()





independent_variables = ["f0_multiplier","jitter_amplitude","jitter_frequency","rd","rpp_k","shimmer_amplitude","shimmer_frequency"]
sound_names = ["female_runn", "female_sust", "male_runn", "male_sust"]
objective_measure = ["cpp", "ps", "pesq", "se", "hnr", "jitter_and_shimmer"]
params_dict = {"f0_multiplier": [0.5, 0.1, 4, 10], "jitter_amplitude": [0, 2, 50, 1], "jitter_frequency": [0, 100, 22000, 1], "rd": [0.35, 0.05, 4, 100], "rpp_k": [0, 0.2, 5, 10], "shimmer_amplitude": [0, 2, 50, 1], "shimmer_frequency": [0, 100, 22000, 1]}
data_leaps = {"female_runn": 16, "female_sust": 24, "male_runn": 13, "male_sust": 28}

for ind_var in independent_variables:
	for sounds in sound_names:
		for measure in objective_measure:
			file_name = ind_var + "_" + measure + "_" + sounds + ".csv"

			param_min = params_dict[ind_var][0]
			param_leap = params_dict[ind_var][1]
			param_max = params_dict[ind_var][2]
			param_multiplier = params_dict[ind_var][3]

			print("Cleaning " + file_name)
			#print("min: " + str(param_min) + " max: " + str(param_max) + " leap: " + str(param_leap))
			if (measure == "cpp"):
				cpp_cleaner("../outputs/analysis/" + file_name,
							"../outputs/cleaned/" + file_name,
							ind_var,
							param_min, param_leap, param_max, param_multiplier)


			if (measure == "ps"):
				ps_cleaner("../outputs/analysis/" + file_name,
							"../outputs/cleaned/" + file_name,
							ind_var,
							param_min, param_leap, param_max, param_multiplier)

			if (measure == "pesq"):
				pesq_cleaner("../outputs/analysis/" + file_name,
							"../outputs/cleaned/" + file_name,
							ind_var,
							param_min, param_leap, param_max, param_multiplier)

			if measure == "se" and (file_name == "female_sust" or file_name == "male_sust") and (ind_var == "rd" or ind_var == "rpp_k" or ind_var == "f0_multiplier"):
				se_cleaner("../outputs/analysis/" + file_name,
							"../outputs/cleaned/" + file_name,
							ind_var,
							param_min, param_leap, param_max, param_multiplier)

			if (measure == "hnr"):
				hnr_jys_cleaner("../outputs/analysis/" + file_name,
							"../outputs/cleaned/" + file_name,
							ind_var,
							param_min, param_leap, param_max, param_multiplier, 1)

			if (measure == "jitter_and_shimmer"):
				hnr_jys_cleaner("../outputs/analysis/" + file_name,
							"../outputs/cleaned/" + file_name,
							ind_var,
							param_min, param_leap, param_max, param_multiplier, data_leaps[sounds])
