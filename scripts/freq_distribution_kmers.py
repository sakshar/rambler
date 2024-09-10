import xlsxwriter
import sys
import numpy as np


def write_histo_to_xcel_get_lower_upper(input, output):
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet()
    
    freqs, no_of_kmers = [], []

    histo = open(input, 'r')
    i = 1
    for line in histo:
        args = line.strip().split()
        worksheet.write("A" + str(i), args[0])
        worksheet.write("B" + str(i), args[1])
        i += 1
        if args[0] == '0':
            continue
        freqs.append(int(args[0]))
        no_of_kmers.append(int(args[1]))

    workbook.close()

    # calculating the boundaries for the frequency ditribution curve

    start_index = 0
    while no_of_kmers[start_index]/no_of_kmers[start_index+1] > 1.25:
        start_index += 1
    # print(f"start index: {start_index}")
    truncated_no_of_kmers = no_of_kmers[start_index:]
    max_index = no_of_kmers.index(max(truncated_no_of_kmers))
    # print(f"max index: {max_index}")
    end_index = start_index + 2*(max_index - start_index)
    # print(f"end index: {end_index}")

    # calculating the mean and the standard deviation of the distribution

    F = np.array(no_of_kmers[start_index:end_index+1])
    X = np.array(freqs[start_index:end_index+1])

    mean = np.sum(X*F)/np.sum(F)
    variance = np.sum(F*((X-mean)*(X-mean)))/np.sum(F)
    stdev = np.sqrt(variance)

    print(f"mean: {mean} std: {stdev}")
    
    # calculating the lower and upper bounds for extracting unikmers
    
    lower = max(int(np.floor(mean - 3*stdev)), 5)
    upper = int(np.ceil(mean + 3*stdev))

    print(f"lower: {lower}")
    print(f"upper: {upper}")
    

if __name__ == "__main__":
    write_histo_to_xcel_get_lower_upper(sys.argv[1], sys.argv[2])
