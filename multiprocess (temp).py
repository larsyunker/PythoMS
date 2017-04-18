#from _classes._Spectrum import Spectrum
from multiprocess import Queue,Process,Pool,cpu_count,freeze_support
freeze_support()

#spec = Spectrum(
#10, # decimal places
#start = None, # no minimum mass
#end = None, # no maximum mass
#empty = True, # whether or not to use emptyspec
#filler = 0., # fill with zeros, not None
#)

#import sys,os
#sys.path.append(os.path.realpath(__file__))

def process_combination(queue):
    while True:
        val = queue.get()
        #spec.addval(val,1.)
        print val




#queue = Queue(cpu_count()) # create Queue instance
queue = Queue()

#pc = Pool(cpu_count(), process_combination,(queue,))
#pc = Process(target=process_combination, args=((queue),)) # inititate the process instance
#pc.start() # start the processing

#procs = []
#for i in range(cpu_count()):
#    pc = Process(target=process_combination, args=((queue),)) # inititate the process instance
#    procs.append(pc)
#    pc.start()
pc = Pool(cpu_count(),process_combination,(queue,))
pc.start()

for i in range(10):
    queue.put(i)


print queue.qsize()
pc.join()

print len(spec)

#pc.join()