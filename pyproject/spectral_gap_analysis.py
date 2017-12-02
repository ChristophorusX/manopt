import numpy as np
import sync_generator as syncgen
import sbm_generator as sbmgen
import burer_monteiro as bm
from sklearn import cluster
import aux


def _spectral_gap(A, z):
    gap = aux.laplacian_eigs(A, z)[1]
    print('Spectral gap: {}'.format(gap))
    return gap


def _gen_sync(n, percentage, snr):
    return syncgen.synchronization_usual(n, percentage, snr)


def _gen_sbm(n, a, b):
    return sbmgen.sbm_logarithm(n, a, b)


def _check_spectral_gap(A, z):
    if _spectral_gap(A, z) > 0:
        return True
    else:
        return False


def search_counter_eg(n, snr, n_iter, n_trail):
    found_target = False
    percentage = .5
    examples = []

    while found_target == False and snr > 3:
        snr -= .1
        print('Starting loops with SNR = {}...'.format(snr))

        for i in range(n_iter):
            print('Loop #{}'.format(i + 1))
            A, z = _gen_sync(n, percentage, snr)
            # A, z = _gen_sbm(n, 10, 2)

            if _check_spectral_gap(A, z):
                print('>>>Found matrix where SDP tight...')
                for j in range(n_trail):
                    print(
                        '>>>>>>Finding global optimizer with BM (trail {})...'.format(j + 1))
                    Q = bm.augmented_lagrangian(
                        A, 2, plotting=False, printing=False)
                    kmeans = cluster.KMeans(
                        n_clusters=2, random_state=0).fit(Q)
                    clustering = 2 * kmeans.labels_ - 1
                    err = aux.error_rate(clustering, z.ravel())
                    print('The error rate for BM is: {}...'.format(err))
                    if err != 0:
                        found_target = True
                        print('One instance found when SNR = {}!'.format(snr))
                        gap = aux.laplacian_eigs(A, z)[1]
                        if gap > .01:
                            example = CounterExample(A, z, Q, gap)
                            examples.append(example)
                            print(A)
    return examples


class CounterExample():

    def __init__(self, A, z, Q, gap):
        self.A = A
        self.z = z
        self.Q = Q
        self.gap = gap

    def get_noise(self):
        return self.A - self.z.dot(self.z.T)


if __name__ == '__main__':
    examples = search_counter_eg(1000, 3.9, 20, 10)
