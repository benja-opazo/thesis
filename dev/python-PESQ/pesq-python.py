from scipy.io import wavfile
from pesq import pesq
import sys, getopt

def main(argv):
	reffile = ''
	degfile = ''
	try:
		opts, args = getopt.getopt(argv,"hr:d:",["ref=","deg="])
	except getopt.GetoptError:
		print('test.py -reg <reffile> -deg <degfile>')
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print('test.py -reg <reffile> -deg <degfile>')
			sys.exit()
		elif opt in ("-r", "--ref"):
			reffile = arg
		elif opt in ("-d", "--deg"):
			degfile = arg

	rate, ref = wavfile.read(reffile)

	rate, deg = wavfile.read(degfile)

	print(pesq(rate, ref, deg,'wb'), end ="")

if __name__ == "__main__":
	main(sys.argv[1:])


