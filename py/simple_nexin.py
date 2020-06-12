import os, sys, math, json, csv

#DEBUG
#inpath = '/media/bcowley/Transcend/flowcar-backup/01/Nexus_data/session_01-01/'
#nxfile = 'nexus_1573639466.0630298.txt'


# Get signal out of JSON file
def getsig(nxfile, sigs):
    nxdata = []
    #    create default data in case there are reading problems
    nx_ts_i = {'ts':math.nan}
    nx_sg_i = dict((k,math.nan) for k in sigs)
    with open(nxfile) as f:
        i = 0
        for line in f:
            i += 1
            try:
                json_object = json.loads(line)
                nx_ts_i = json_object[0]
                nx_sg_i = json_object[1]
            except:
                print("JSON line-" + str(i) + " read error for line: " + line)
            else:
                fields = [nx_ts_i[1]]
                for s in sigs:
                    fields.append(nx_sg_i[s])
                nxdata.append(fields)
            
        f.close()
    return nxdata


# Function walks through folders under given path & parses any nexus files found
# default-coded name of signal to read
def parse_nxs(inpath, sigs):
    sigmap = {'EDA':'E', 'EOG':['A', 'B'], 'BVP':'F'}
    inpath = os.path.abspath(inpath)
    for root, subdirs, files in os.walk(inpath):
        for f in files:
            if f.startswith('nexus'):
                nxf = os.path.join(root, f)
                print('Found Nexus, parsing ' + sigs + ' from ' + nxf)
                dat = getsig(nxf, sigmap[sigs])
                outf = os.path.join(root, sigs + '_' + f[0:-4] + '.csv')
                outf = open(outf, "w", newline="")
                writer = csv.writer(outf)
                writer.writerows(dat)

            if subdirs != []:
                print('Going deeper')
                for s in subdirs:
                    parse_nxs(s)

if __name__ == "__main__":
    parse_nxs(sys.argv[1], sys.argv[2])

#DEBUG
#parse_nxs(inpath, 'EOG')
