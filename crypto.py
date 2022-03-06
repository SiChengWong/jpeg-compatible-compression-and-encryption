from struct import unpack
import rc4

from decoder import Stream
from decoder import GetArray
from decoder import HuffmanTable
from decoder import RemoveFF00

class CryptoPretreatment:

    def byteToBitList(self, byte: int) -> list[int]:
        bit_list = []
        for i in range(8):
            bit_list.insert(0, byte & 1)
            byte = byte >> 1
        return bit_list

    def bitListToByte(self, bit_list: list[int]) -> int:
        byte = 0
        for i in range(len(bit_list)):
            byte = (byte << 1) + bit_list[i]
        for i in range(len(bit_list), 8):
            byte = byte << 1
        return byte

    def bitsLenSort(self, amplitudes: list[list[int]]) -> list[list[int]]:
        """
        sort the array 'amplitudes' in method 'EncryptMatrix' of class 'Crypto'
        in descending order of bit length(namely, code value)
        :param amplitudes: the array
        :return: sorted array
        """
        bits_len_dict = {}
        for i in range(len(amplitudes)):
            if amplitudes[i][1] not in bits_len_dict:
                bits_len_dict[amplitudes[i][1]] = []
            bits_len_dict[amplitudes[i][1]].append(amplitudes[i])
        sorted_amplitudes = []
        for i in sorted(bits_len_dict, reverse=True):
            sorted_amplitudes += bits_len_dict[i]
        return sorted_amplitudes

    def divideAmplitudesIntoBytes(self, amplitudes: list[list[int]]) -> list[int]:
        """
        append bits with variable length into a list of integer ranging from 0 to 255,
        in order to operate RC4 encrypting conveniently
        :param amplitudes: list containing bits information,
                           defined in method 'EncryptMatrix' of class 'Crypto'
        :return: a list of integer ranging from 0 to 255
        """
        bit_list = []
        for i in range(len(amplitudes)):
            bits = amplitudes[i][2]
            bit_list_of_one_elem = []
            while len(bit_list_of_one_elem) < amplitudes[i][1]:
                bit_list_of_one_elem.insert(0, bits & 1)
                bits = bits >> 1
            bit_list += bit_list_of_one_elem
        # padding for dividing into bytes
        while len(bit_list) % 8 != 0:
            bit_list.append(0)
        byte_list = []
        for i in range(len(bit_list) >> 3):
            byte_list.append(self.bitListToByte(bit_list[i << 3: (i<<3) + 8]))
        return byte_list

    def modifyBitInStream(self, st: Stream, pos: int, bit: int) -> None:
        """
        modify the 'pos'-th bit in 'st' according to 'bit'
        :param st: target Stream object
        :param pos: position of bit to be modified
        :param bit:
        :return:
        """
        byte = st.data[pos >> 3]
        # divide byte into bits
        bit_list = []
        for i in range(8):
            bit_list.insert(0, byte & 1)
            byte = byte >> 1
        bit_list[pos & 7] = bit
        st.data[pos >> 3] = 0
        for i in range(8):
            st.data[pos >> 3] = (st.data[pos >> 3] << 1) + bit_list[i]


    def modifyAmplitudes(self, st: Stream, amplitudes: list[list[int]], encrypted_byte_list: list[int]) -> None:
        """
        use encrypted_bit_list to rewrite amplitudes
        """
        encrypted_bit_list = []
        for i in range(len(encrypted_byte_list)):
            encrypted_bit_list += self.byteToBitList(encrypted_byte_list[i])
        k = 0   # index of 'encrypt_bit_list'
        for i in range(len(amplitudes)):
            for j in range(amplitudes[i][1]):
                self.modifyBitInStream(
                    st,
                    amplitudes[i][0] + j,
                    encrypted_bit_list[k]
                )
                k += 1

class Crypto:
    """
    open a JPEG image and encrypt it
    """

    def __init__(self, image_file: str, K: list):
        """
        set value to necessary parameters
        :param image_file: name of image file to be encrypted
        :type K: list of integer range from 0 to 255
        :param K: secrete key to encrypt the image
        """
        self.K = K
        self.rc4 = rc4.RC4(self.K)
        self.huffman_tables = {}
        self.image_file = image_file
        with open(image_file, "rb") as f:
            self.img_data = f.read()

    def BaselineDCT(self, data):
        """
        get basic information of JPEG image file, like height and width
        :param data: bit stream
        :return: None
        """
        hdr, self.height, self.width, components = unpack(">BHHB", data[0:6])
        print("size %ix%i" % (self.width, self.height))

    def decodeHuffman(self, data):
        while len(data) > 0:
            offset = 0
            # bit 0..3: number of HT (0..3, otherwise error)
            # bit 4: type of HT, 0 = DC table, 1 = AC table
            (header,) = unpack("B", data[offset: offset + 1])
            print(header, header & 0x0F, (header >> 4) & 0x0F)
            offset += 1

            lengths = GetArray("B", data[offset: offset + 16], 16)
            offset += 16

            elements = []
            for i in lengths:
                elements += GetArray("B", data[offset: offset + i], i)
                offset += i

            hf = HuffmanTable()
            hf.GetHuffmanBits(lengths, elements)
            self.huffman_tables[header] = hf
            data = data[offset:]

    def EncryptMatrix(self, st: Stream, idx: int) -> None:
        # array 'amplitudes' record start position of amplitudes
        # of DC coefficients and non-zero AC coefficients,
        # in which each element is a tuple (start_position, bits),
        # whose:
        # 1st component denotes the amplitude start position in st,
        # 2nd one denotes the amplitude length in bits,
        # 3rd one denotes bit stream representing amplitude value
        amplitudes = []

        # get DC coefficient
        code = self.huffman_tables[0 + idx].GetCode(st) # SIZE
        bits = st.GetBitN(code)                         # raw AMPLITUDE bits
        amplitudes.append([st.pos, code, bits])

        l = 1
        while l < 64:
            code = self.huffman_tables[16 + idx].GetCode(st)
            if code == 0:
                # End of Block
                break
            if code > 15:
                # 4 MSB denotes zero run-length
                l += code >> 4
                # 4 LSB denotes amplitude length in bits
                code = code & 0x0F
            bits = st.GetBitN(code)
            amplitudes.append([st.pos, code, bits])
            if l < 64:
                l += 1

        # sort AC coefficient amplitudes in descending order of code
        # DC coefficient amplitude remains still
        pretreatment = CryptoPretreatment()
        amplitudes = [amplitudes[0]] + pretreatment.bitsLenSort(amplitudes[1:])
        # divide bit stream into list of byte
        byte_list = pretreatment.divideAmplitudesIntoBytes(amplitudes)
        # encrypt list of byte with RC4
        encrypted_byte_list = self.rc4.encrypt(byte_list)
        pretreatment.modifyAmplitudes(st, amplitudes, encrypted_byte_list)

    def StartOfScan(self, data, hdrlen):
        data, lenchunk = RemoveFF00(data[hdrlen:])

        # data was assigned directly to st.data,
        # so data can be modified by set values in st.data
        st = Stream(data)
        for y in range(self.height // 8):
            for x in range(self.width // 8):
                self.EncryptMatrix(st, 0)
                self.EncryptMatrix(st, 1)
                self.EncryptMatrix(st, 1)


        return

    def encrypt(self):
        """
        encrypt JPEG image and save as a new file
        :return: cipher text
        """
        data = self.img_data
        cipher_bytes = bytes(0)
        while True:
            (marker,) = unpack(">H", data[0:2])
            # add marker in clear
            cipher_bytes += data[0:2]
            if marker == 0xFFD8:
                data = data[2:]
            elif marker == 0xFFD9:
                return
            else:
                (len_chunk,) = unpack(">H", data[2:4])
                len_chunk += 2
                chunk = data[4:len_chunk]
                if marker == 0xFFC4:
                    self.decodeHuffman(chunk)
                    # add chunk length and content
                    cipher_bytes += data[2:len_chunk]
                elif marker == 0xFFC0:
                    self.BaselineDCT(chunk)
                    # add chunk length and content
                    cipher_bytes += data[2:len_chunk]
                elif marker == 0xFFDA:
                    len_chunk = self.StartOfScan(data, len_chunk)
                else:
                    cipher_bytes += data[2:len_chunk]
                data = data[len_chunk:]
            if len(data) == 0:
                break