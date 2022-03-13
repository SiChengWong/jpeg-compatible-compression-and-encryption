import base64


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
        self.initialize(self.K)

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

    def encrypt(self, M: list, continue_flag = False):
        """
        encrypt message M with RC4 algorithm
        :param continue_flag: if continue the encryption without re-initializing the S-box
        :type M: list of byte
        :param M: plain text
        :return: cipher text
        """
        if continue_flag == False:
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

    def decrypt(self, C: list, continue_flag = False):
        """
        encrypt message M with RC4 algorithm
        :param continue_flag: if continue the encryption without re-initializing the S-box
        :type C: list of byte
        :param C: cipher text
        :return: plain text
        """
        if continue_flag == False:
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


def intToBitList(n: int) -> list[int]:
    bit_list = []
    while n > 0:
        bit_list.insert(0, n & 1)
        n = n >> 1
    return bit_list


def paddingToLen(l: list[int], length: int, pos: int = 0):
    while len(l) < length:
        l.insert(pos, 0)
    return l


class RC4RandBitGenerator:
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

    def initialize(self, key: list[int]):
        self.S = self.T = []
        for i in range(256):
            self.S.append(i)
            self.T.append(key[i % len(key)])
        j = 0
        for i in range(256):
            j = (j + self.S[i] + self.T[i]) % 256
            self.swap(self.S, i, j)

    def __init__(self, key: list[int]):
        self.key = [key[i] & 0xFF for i in range(len(key))]
        self.initialize(self.key)

        self.i = self.j = 0
        self.bit_list = []

    def genRandBit(self) -> int:
        if len(self.bit_list) == 0:
            self.i = (self.i + 1) & 0xFF                        # increase by 1, and then mod 256
            self.j = (self.j + self.S[self.i]) & 0xFF           # increase by self.S[self.i], and then mod 256
            self.swap(self.S, self.i, self.j)
            t = (self.S[self.i] + self.S[self.j]) & 0xFF
            self.bit_list += paddingToLen(intToBitList(self.S[t]), 8)
        return self.bit_list.pop(0)

def testRC4():
    M = r"HelloWorld"
    K = r"SecreteKey"

    rc4 = RC4([ord(k) for k in K])

    # encrypt message and display
    C = rc4.encrypt([ord(m) for m in M])
    cipher_text = base64.encodebytes(bytes(C))
    print(cipher_text)

    # decrypt cipher text and display
    M_ = rc4.decrypt(C)
    print("".join([chr(m_) for m_ in M_]))

def testRC4randGenerator():
    key = [1,2,3,4,5]
    rc4_rand_bit_generator = RC4RandBitGenerator(key)
    for i in range(1000):
        print(rc4_rand_bit_generator.genRandBit())

if __name__ == '__main__':
    testRC4randGenerator()