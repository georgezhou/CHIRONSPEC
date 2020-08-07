import os,sys

def main(user):

    cwd = os.getcwd()
    os.chdir("/data/tfop/chiron_data/")

    cmd = "wget -nH --cut-dirs=1 --recursive -np -R 'index.html' -nc --user=smarts --password='02Babg*!' ftp://ftp.chara.gsu.edu:2997/"+user+"/"


    os.system(cmd)
    os.chdir(cwd)


if __name__ == "__main__":
    main("zhou")

    
