import time

def start_time(string):
	print(string)
	t0 = time.time()
	return t0

def end_time(t0):
	print('Done. (%.2f s)\n' % (time.time() - t0))