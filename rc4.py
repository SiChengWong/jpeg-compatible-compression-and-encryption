class RC4:
    """
    RC4 encryption and decryption algorithm
    """

    def __init__(self, K: list):
        """
        initialization
        :type K: list of integer range from 0-255
        :param K: secrete key
        """
        self.K = K

    def swap(self, array: list, i: int, j: int):
        """
        swap elements with indices i,j in a list
        :param array: list to be permuted
        :param i: one of the indices to be swapped
        :param j: another index to be swapped
        :return: None
        """
        temp = array[i]
        array[i] = array[j]
        array[j] = temp

    def initialize(self, K):
        self.S = self.T = []
        for i in range(256):
            self.S.append(i)
            self.T.append(K[i % len(K)])
        j = 0
        for i in range(256):
            j = (j + self.S[i] + self.T[i]) % 256
            self.swap(self.S, i, j)

    def encrypt(self, M: list):
        """
        encrypt message M with RC4 algorithm
        :type M: list of byte
        :param M: plain text
        :return: cipher text
        """
        self.initialize(self.K)
        i = j = 0
        C = []
        for m in M:
            i = (i + 1) % 256
            j = (j + self.S[i]) % 256
            self.swap(self.S, i, j)
            t = (self.S[i] + self.S[j]) % 256
            C.append(m ^ self.S[t])
        return C

    def decrypt(self, C: list):
        """
        encrypt message M with RC4 algorithm
        :type C: list of byte
        :param C: cipher text
        :return: plain text
        """
        self.initialize(self.K)
        i = j = 0
        M = []
        for c in C:
            i = (i + 1) % 256
            j = (j + self.S[i]) % 256
            self.swap(self.S, i, j)
            t = (self.S[i] + self.S[j]) % 256
            M.append(c ^ self.S[t])
        return M

def main():
    M = r"HelloWorld"
    K = r"SecreteKey"

    rc4 = RC4([ord(k) for k in K])

    # encrypt message and display
    C = rc4.encrypt([ord(m) for m in M])
    print("".join([chr(c) for c in C]))

    # decrypt cipher text and display
    M_ = rc4.decrypt(C)
    print("".join([chr(m_) for  m_ in M_]))

if __name__ == '__main__':
    main()