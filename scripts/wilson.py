from math import sqrt
def wilson(k, n, z=1.959963985):
    if n == 0: return None
    p = k/n
    d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = (z/d)*sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (100*max(0.0, c-h), 100*min(1.0, c+h))

def fmt(label, k, n):
    lo, hi = wilson(k, n)
    print(f"{label:<34} {k}/{n} = {100*k/n:6.2f}%  CI {lo:5.1f}-{hi:5.1f}")

print("== DISCOVERY (602 labeled; 19 L, 583 non-L; 11 FP) ==")
fmt("sensitivity", 19, 19)
fmt("specificity", 583-11, 583)
fmt("PPV", 19, 30)
fmt("NPV", 572, 572)

print("\n== EXTERNAL n=301, gold standard ST6 (22 L, 279 non-L, 2 FP) ==")
fmt("sensitivity", 22, 22)
fmt("specificity", 279-2, 279)
fmt("PPV", 22, 24)

print("\n== EXTERNAL n=301, gold standard <=50 alleles (23 L, 278 non-L, 1 FP) ==")
fmt("sensitivity", 23, 23)
fmt("specificity", 278-1, 278)
fmt("PPV", 23, 24)

print("\n== EXTERNAL n=295, ST6 (22 L, 273 non-L, 2 FP) ==")
fmt("sensitivity", 22, 22)
fmt("specificity", 273-2, 273)
fmt("PPV", 22, 24)

print("\n== EXTERNAL n=295, <=50 alleles (23 L, 272 non-L, 1 FP) ==")
fmt("sensitivity", 23, 23)
fmt("specificity", 272-1, 272)
fmt("PPV", 23, 24)

print("\n== ST3 marker-carrying proportions, Fisher exact ==")
try:
    from scipy.stats import fisher_exact
    table = [[12, 38-12], [1, 27-1]]
    odds, p = fisher_exact(table, alternative='two-sided')
    print(f"discovery 12/38 = {100*12/38:.1f}%   external 1/27 = {100*1/27:.1f}%")
    print(f"Fisher two-sided p = {p:.4f}   odds ratio = {odds:.2f}")
except ImportError:
    print("scipy missing")
