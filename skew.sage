R.<q,t> = QQ['q,t']
Sym = SymmetricFunctions(FractionField(R))
P = Sym.macdonald().P()
Q = Sym.macdonald().Q()

from sage.combinat.q_analogues import q_pochhammer

def struc_cons(l,m,n):
	return (P(m)*P(n)).scalar_qt(Q(l))

def nl(l):
	return sum(i*l[i] for i in range(len(l)))

def nlm(l,m):
	lp = Partition(l).conjugate()
	mp = Partition(m).conjugate()
	return sum(binomial(lp.get_part(i)-mp.get_part(i),2) for i in range(len(lp)))

def bl(l):
	l = Partition(l)
	return prod((1-q^(l.arm_length(i,j)) * t^(l.leg_length(i,j)+1)) / (1- q^(l.arm_length(i,j)+1) * t^(l.leg_length(i,j))) for i,j in l.cells())

def Pprincipal(l):
	l = Partition(l)
	return t^(nl(l)) * prod(1/(1-q^(l.arm_length(i,j)) * t^(l.leg_length(i,j)+1)) for i,j in l.cells())

def Qlm(l,m):
	l = Partition(l)
	m = Partition(m)
	return sum(struc_cons(l,m,nu) * bl(nu) * Pprincipal(nu) for nu in Partitions(sum(l)-sum(m), outer=l))

def Qratio(l,m):
	return Qlm(l,m)/(bl(l)*Pprincipal(l))

def Plm(l,m):
	return bl(m) * Qlm(l,m) / bl(l)

def Pratio(l,m):
	return bl(m) * Qratio(l,m)

def QHLprincipal(l,m):
	q,t = var('q,t')
	l = Partition(l)
	m = Partition(m)
	lp = l.conjugate()
	mp = m.conjugate()
	d = m.to_exp()
	return t^nlm(l,m) * prod(q_pochhammer(d[i], t^(1+lp.get_part(i)-mp.get_part(i)), t)/q_pochhammer(d[i], t, t) for i in range(m[0]))