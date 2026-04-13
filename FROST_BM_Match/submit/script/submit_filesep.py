#~/usr/bin/env python3

import subprocess

filepath = '../shell/run_filesep.sh'

nmonth = 4

def main():
    for i in range (0, nmonth) :
        year, month, startday, endday = get_input_date(i)
        for day in range (startday, endday + 1) :
#           print (year, month, day)
#           subprocess.call(["bsub -q s " + filepath + " " + str(year) + " " + str(month) + " " + str(day)], shell=True)
            subprocess.call([filepath + " " + str(year) + " " + str(month) + " " + str(day)], shell=True)

def get_input_date(i):
    if i == 0 :
        year = 2019
        month = 11
        startday = 7
        endday = 30
    elif i == 1 :
        year = 2019
        month = 12
        startday = 1
        endday = 19
    elif i == 2 :
        year = 2020
        month = 1
        startday = 14
        endday = 31
    elif i == 3 :
        year = 2020
        month = 2
        startday = 1
        endday = 12

    return year, month, startday, endday
            
if __name__ == "__main__":
    main()
