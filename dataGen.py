import numpy as np

def dataGen(file_name, bound = 1, dim = 2, n_cluster = 5, n_samples = 20000):
    print(dim, n_cluster, n_samples)
    f = open(file_name,'w')
    uniform_ratio = 0.3
    f.write(str(int(n_samples*(1.+uniform_ratio)))+" "+str(dim)+"\n")
    for x in range(int(n_samples*uniform_ratio)):
        s = np.random.uniform(-2*bound, 2*bound, dim)
        line = ""
        for i in range(len(s)):
            line += str(s[i])+ " "
        f.write(line+"\n")
    p_bound = 2*bound

    # range = -bound ~ bound
    for x in range(n_samples):
        if x % (n_samples/n_cluster) ==0:
            mu = np.random.uniform(-bound, bound, dim)
            sigma = np.random.uniform(0, bound, dim)*0.1 +0.05
        s = np.random.normal(mu, sigma)
        line = ""
        for i in range(len(s)):
            if s[i] < -p_bound:	s[i] = -p_bound
            if s[i] > p_bound: s[i] = p_bound
            line += str(s[i])+ " "
        f.write(line+"\n")
    f.close()
    return

if __name__ == '__main__':
    bound = 1
    dim = 2
    n_cluster = 5
    n_samples = 200000
    file_name = 'synthetic_data_%d_%d.txt' % (n_cluster, n_samples)
    dataGen(file_name, bound, dim, n_cluster, n_samples)
