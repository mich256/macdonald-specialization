R.<t> = QQ['t']

def skew_schur_ratio(l,mu):
	temp = 0
	l = Partition(l)
	mu = Partition(mu)
	n = len(l)
	for tab in SemistandardTableaux(mu, max_entry=n):
		temp += prod(prod(t^(i-n)*(1-t^(l[n-i]-k+j)) for j,k in tab.cells_containing(i)) for i in range(1,n+1))
	return temp

def nl(l):
	return sum(i*l[i] for i in range(len(l)))

def schur_principal(l):
	l = Partition(l)
	return t^(nl(l))/prod(1-t^(l.hook_length(i,j)) for i,j in l.cells())