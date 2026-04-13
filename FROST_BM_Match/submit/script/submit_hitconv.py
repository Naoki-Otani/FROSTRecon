from submit_filesep import get_input_date
import subprocess

filepath = '../shell/run_hitconv.sh'

nmonth = 4

def main() :
    for i in range (0, nmonth) :
        year, month, startday, endday = get_input_date(i)
        for day in range (startday, endday + 1) :
#           print (year, month, day)
#           subprocess.call(["bsub -q s" + filepath + " " + str(year) + " " + str(month) + " " + str(day)], shell=True)
            subprocess.call([filepath + " " + str(year) + " " + str(month) + " " + str(day)], shell=True)

if __name__ == "__main__" :
    main()
