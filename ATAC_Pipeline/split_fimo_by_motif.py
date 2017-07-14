import sys

fimo = open(str(sys.argv[1]))

# dataToDictionary-------------------------------------------------------------
#
# Represent fimo file as a dictionary in which the key is the motif name and
# the value is the rest of the line.
#------------------------------------------------------------------------------
def  dataToDictionary(myFile):

	fimo_dict = {}

	file_index = 0
	for line in myFile:
		#if this isn't the first line:
		if file_index > 0:
			#split the line by tabs:
			line_list = line.rsplit('\t')

			#dictionary key is the motif name
			key = line_list[0]

			#split the chromosome number, start, and stop locations
			chr_start_stop = line_list[1].rsplit(':')
			chrm = chr_start_stop[0]
			start_stop = chr_start_stop[1].rsplit('-')

			#compute actual start and stop:
			start = str(int(start_stop[0]) + int(line_list[2]))
			stop = str(int(start_stop[0])+ int(line_list[3]))

			#create dictionary entry
			item = [chrm, start, stop, key, line_list[5], line_list[4], line_list[6], line_list[len(line_list)-1][:len(line_list[len(line_list)-1])-1]]
			if key not in fimo_dict.keys():
				fimo_dict[key] = [item]
			else:
				fimo_dict[key].append(item)
		file_index += 1

	return fimo_dict

# splitFimo--------------------------------------------------------------------
#
# split output from fimo into separate files based on their motif name. This is
# done using the helper function dataToDictionary, which creates a dictionary
# from fimo's output in which the keys are the motif names.
#------------------------------------------------------------------------------
def splitFimo():

	#call dataToDictionary:
	my_dict = dataToDictionary(fimo)

	#for each motif:
	for key in my_dict.keys():

		my_items = my_dict[key]
		filename = key + ".motif.bed"
		output_file = open(filename, 'w')

		#write all corresponding items to its output file:
		for item in my_items:
			out_string = ""
			for i in range(len(item)):
				if i < len(item)-1:
					out_string += item[i] + '\t'
				else:
					out_string += item[i] + '\n'
			output_file.write(out_string)

splitFimo()