# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import xlrd

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


def main():
    wb = xlrd.open_workbook('/Users/liyaqi/Documents/生信/Inhouse_cohorts_genes_Version_8_MRR_诊断.xlsx')
    sh = wb.sheet_by_name('Sheet1')
    print(sh.nrows)
    print(sh.ncols)
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
