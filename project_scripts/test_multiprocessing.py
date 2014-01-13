'''
example of mutliprocessing module - add to values contained in a tuple 

'''

import multiprocessing
import cmap.util.progress as progress 

def _add_two(params):
    '''
    input: the index of a signature from the input gctx, output: saves a self connection graph
    '''
    val1 = params[0]
    val2 = params[1]
    sum1 = val1 + val2 
    time.sleep(10)
    return sum1


# make list of tuples
tupList = []
for x1 in np.arange(3):
    for x2 in np.arange(3):
        tup1 = (x1,x2)
        tupList.append(tup1) #add to list of tuples 

# instantiate a progress object
prog = progress.DeterminateProgressBar('self-connection graph builder')
#build graphs in parallel
n_procs=4
pool = multiprocessing.Pool(n_procs)
rs = pool.map_async(_add_two,tupList)
pool.close() # No more work
while (True):
    if (rs.ready()): break
    remaining = rs._number_left
    prog.show_message('Waiting for {0} tasks to complete...'.format(remaining))
    time.sleep(0.1)
rs.get()

### attempt 2
import sys, time, random, multiprocessing

def talk(n1,n2):
    for i in xrange(n1,n2):
        print "hello from %s" % i
        time.sleep(random.random())  #pretend to do some work

p1 = multiprocessing.Process(target=talk, args=(2,5,))
p2 = multiprocessing.Process(target=talk, args=(7,9,))
p1.start(); p2.start()
p1.join() ; p2.join()    


# pool = multiprocessing.Pool()
# tupList = [(2,5),(7,9)]
# rs = pool.map_async(talk,tupList)


# ### attempt 3
# from multiprocessing import Pool

# def f(x):
#     return x*x


# if __name__ == '__main__':
#     pool = Pool(processes=4)              # start 4 worker processes
#     result = pool.apply_async(f, [10])    # evaluate "f(10)" asynchronously
#     print result.get(timeout=1)           # prints "100" unless your computer is *very* slow
#     print pool.map(f, range(10))          # prints "[0, 1, 4,..., 81]"


# pool2 = Pool(processes=4)
# result = pool.map_async(talk,[(2,5,)])
# result = pool.map_async(f,[4,9])


# ### attempt 4 - use corey's q_runner

import cmap.util.queue as cq

# def _worker1(input, output):
#     '''example worker process'''
#     for c in iter(input.get, 'STOP'):
#         # do some stuff
#         subprocess.Popen('sleep 1', shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
#         # stick the ouput somewhere
#         output.put(str(c) + 'a')

# res1 = cq.q_runner(range(10), _worker1, n_procs=2)

# res1 = cq.q_runner([(2),(3)], f, n_procs=2)


o1 = _svm_worker(zFrm, brd, probeIDs)
res1 = cq.q_runner(tupList, _svm_worker, n_procs=6)
    


def _worker1(input, output):
    '''example worker process'''
    for c in iter(input.get, 'STOP'):
        # do some stuff
        subprocess.Popen('sleep 1', shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        # stick the ouput somewhere
        output.put(str(c) + 'a')
 
 
if __name__ == '__main__':
    n_procs = 6
    list_item = range(10)
    r1 = cq.q_runner(list_item, _worker1, n_procs=n_procs)
    print 'first process completed'
    r2 = cq.q_runner(list_item, _worker1, n_procs=n_procs)
    print 'second process completed'  
    print r2    

# ###



