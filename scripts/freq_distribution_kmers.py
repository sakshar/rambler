import xlsxwriter
import sys

def write_histo_to_xcel(input, output):
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet()

    histo = open(input, 'r')
    i = 1
    for line in histo:
        args = line.strip().split()
        worksheet.write("A" + str(i), args[0])
        worksheet.write("B" + str(i), args[1])
        i += 1
    workbook.close()


if __name__ == "__main__":
    write_histo_to_xcel(sys.argv[1], sys.argv[2])