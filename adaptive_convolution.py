import numpy as np
import math as m


class AdaptiveConvolution:
    def __init__(self, map, k, exp, l_0, sigma_0, method='round_gauss'):
        # k is kernel size, must be odd, n is exponent of the adaptive factor
        self.kernel_function = self.__make_kernel_function(method, l_0, sigma_0, exp)
        self.map = map
        self.k = k
        self.n = map.shape[0]
        self.m = int(self.n+k-1)
        # print(self.n, self.k, self.m)
        self.padded_map = np.zeros((self.m, self.m))
        for i in range(self.m):
            for j in range(self.m):
                if (self.k - 1)/2 - 1 < i < self.n + (self.k - 1)/2 and \
                        (self.k - 1) / 2 - 1 < j < self.n + (self.k - 1) / 2:
                    self.padded_map[j][i] = map[int(j-(self.k-1)/2)][int(i-(self.k-1)/2)]
        self.conv_map = np.zeros((self.n, self.n))

    def __make_kernel_function(self, method, l_0, sigma_0, n):
        if method == 'round_gauss':
            if n == 0:
                def _kernel_function(val, x, y):
                    l = l_0
                    return m.exp(-(x*x+y*y)/(l*l))
                return _kernel_function
            else:
                def _kernel_function(val, x, y):
                    l = l_0 * m.pow(m.fabs(val)/sigma_0, n)
                    return m.exp(-(x*x+y*y)/(l*l))
                return _kernel_function
        elif method == 'round_exp':
            if n == 0:
                def _kernel_function(val, x, y):
                    l = l_0
                    return m.exp(-m.sqrt((x*x+y*y))/l)
                return _kernel_function
            else:
                def _kernel_function(val, x, y):
                    l = l_0 * m.pow(m.fabs(val) / sigma_0, n)
                    return m.exp(-m.sqrt((x*x+y*y))/l)
                return _kernel_function

    def convolve(self):
        for i in range(self.n):
            for j in range(self.n):
                self.make_kernel(self.map[j][i])
                for u in range(self.k):
                    for v in range(self.k):
                        self.conv_map[j][i] += self.padded_map[j+v][i+u] * self.kernel[v][u]

    def make_kernel(self, val):
        k = self.k
        self.kernel = np.zeros((k, k))
        for i in range(- int((k - 1) / 2), int((k - 1) / 2 + 1)):
            for j in range(- int((k - 1) / 2), int((k - 1) / 2 + 1)):
                self.kernel[int(j+(k-1)/2)][int(i+(k-1)/2)] = self.kernel_function(val, i, j)
        tmp = 0
        for i in range(- int((k - 1 ) / 2), int((k - 1) / 2 + 1)):
            for j in range(- int((k - 1) / 2), int((k - 1) / 2 + 1)):
                tmp += self.kernel[int(j+(k-1)/2)][int(i+(k-1)/2)]
        self.kernel /= tmp


# map = np.ones((10,10))
# print(map)
# conv = AdaptiveConvolution('round_gauss', map, 7, 0, 2, 1e-3 )
# conv.make_kernel(1)
# tmp = 0
# for i in range(conv.kernel.shape[0]):
#    for j in range(conv.kernel.shape[0]):
#        tmp += conv.kernel[j][i]
# print(tmp)
# conv.convolve()
# np.set_printoptions(precision=3)
# print(conv.padded_map)
# print(conv.conv_map)
