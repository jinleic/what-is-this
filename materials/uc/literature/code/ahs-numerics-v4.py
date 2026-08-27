import numpy as np

PHI = (5**0.5-1)/2
I1_START = .60
I1_END = .77
I1_STEP = 1/200
I2_START = .76
I2_END = .98
I2_STEP = 1e-4
I2_MARGIN = 2e-3


def L(x):
	term1 = 2*PHI*(1-x*x)*np.log(1/(x*x)-1)
	term2 = -4*PHI
	term3 = -2*x*x*np.log(x)
	term4 = 2*(x*x-1)*np.log(1-x)
	term5 = x + 2*np.log(x) + 1
	return term1 + term2 + term3 + term4 + term5

def L_table():
	pts = list(np.arange(I1_START,I1_END,I1_STEP))
	print("x \t L(x)")
	for pt in pts:
		print("%.3f \t %.4f" % (pt, L(pt)))

def H(x):
	return -x * np.log(x) - (1-x) * np.log(1-x)
def g1(x):
	return PHI * H(x*x)
def g2(x):
	return x * H(x)

def find_next_pt(pt1):
	# input: pt1, a multiple of I2_STEP
	# output: pt2, the smallest multiple of I2_STEP near pt1 such that g1(pt1) - g2(pt2) > I2_MARGIN
	g1pt1 = g1(pt1)
	pt2 = pt1
	while g2(pt2) < g1pt1 - I2_MARGIN:
		pt2 -= I2_STEP
	return pt2 + I2_STEP

def g_table():
	pts = [I2_END]
	while pts[-1] > I2_START:
		pts.append(find_next_pt(pts[-1]))
	pts.reverse()
	print("x \t g1(x) \t g2(x)")
	for pt in pts:
		print("%.4f \t %.4f \t %.4f" % (pt, g1(pt), g2(pt)))

L_table()
g_table()
