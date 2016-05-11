import sys
import os.path
import getopt
import numpy as np

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
             raise Usage(msg)
        # more code, unchanged
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())

n_cpu=2
n_core=6

input_dir='/Volumes/Q6/NGFS/DECam/stacks'
input_tile=[1,2,3,4,5,6,7,10,13]
input_filter=['u','g','i']

im_file = [input_dir+'/ss_fornax_tile%s_%s_long.003.fits' % (tile,filter) for tile in input_tile for filter in input_filter]
cat_file = [item.replace('.fits','.ldac') for item in im_file]
seg_im_file = [item.replace('.fits','.SEGMENTATION.fits') for item in im_file]

n_im_file=im_file.size

for item in im_file:
	if os.path.isfile(item):
		print 'Processing image '+item

    n_processes=n_cpu*n_core
    n_processes=np.minimum(n_processes,n_im_file)

    tic=time.time()
    j_step=np.int(np.ceil( n_im_file*1./n_processes ))
    j_range=range(0,n_im_file,j_step)
    j_range.append(n_im_file)

    work_queue = multiprocessing.Queue()
    for j in range(np.size(j_range)-1):
      if work_queue.full():
        print "Oh no! Queue is full after only %d iterations" % j
      work_queue.put( (im_file, cat_file, seg_im_file) )

    # create a queue to pass to workers to store the results
    result_queue = multiprocessing.Queue()
    procs=[]

    # spawn workers
    for j in range(n_processes):
      worker = Worker_sex(work_queue, result_queue)
      procs.append(worker)
      worker.start()
      time.sleep(30)

    for j in range(n_processes):
      result_queue.get()

    for p in procs:
      p.join()


"""
f = open('survey_calib.dat', 'r')  # We need to re-open the file

data = []
for line in f:
  line = line.strip()
  columns = line.split()
  if columns[0] != '#':
    source = {}
    source['im_file'] = columns[0]
    source['mjd'] = float(columns[1])
    source['object'] = columns[2]
    source['filter'] = columns[3]
    data.append(source)

f.close()
"""
