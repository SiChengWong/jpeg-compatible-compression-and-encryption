from struct import unpack
from struct import pack
from rc4 import RC4


def remove0xFF00(l: list[int]) -> list[int]:
    """
    Removes 0x00 after 0xff in the image scan section of JPEG
    :param l: list to be removed 0xFF00
    :return: list after removing 0xFF00
    """
    l_without_0xFF00 = [l[0]]
    for i in range(1, len(l)):
        if l[i] == 0x00:
            if l[i - 1] != 0xFF:
                l_without_0xFF00.append(l[i])
        else:
            l_without_0xFF00.append(l[i])
    return l_without_0xFF00


def add0xFF00(l: list[int]) -> list[int]:
    l_added_0xFF00 = []
    for i in range(len(l)):
        l_added_0xFF00.append(l[i])
        if l[i] == 0xFF:
            l_added_0xFF00.append(0x00)
    return l_added_0xFF00


def GetArray(elem_type, l, length):
    """
    A convenience function for unpacking an array from bitstream
    """
    s = ""
    for i in range(length):
        s = s + elem_type
    return list(unpack(s, l[:length]))


def byteListToInt(l: list[int]) -> int:
    """
    convert list of integer ranging from 0 to 255 into integer
    :param l: list
    :return: result
    """
    res = 0
    byte_list = bytes(l)
    for i in range(len(byte_list)):
        res = (res << 8) + byte_list[i]
    return res


def intToBitList(n: int) -> list[int]:
    bit_list = []
    while n > 0:
        bit_list.insert(0, n & 1)
        n = n >> 1
    return bit_list


def bitListToInt(bit_list: list[int]) -> int:
    n = 0
    for b in bit_list:
        n = (n << 1) + b
    return n


def bitListToByteList(bit_list: list[int]) -> list[int]:
    while len(bit_list) % 8 != 0:
        bit_list.append(0)
    byte_list = []
    for i in range(len(bit_list) // 8):
        byte_list.append(bitListToInt(bit_list[i * 8: i * 8 + 8]))
    return byte_list


def byteListToBitList(byte_list: list[int]) -> list[int]:
    bit_list = []
    for byte in byte_list:
        bit_list += paddingToLen(intToBitList(byte), 8)
    return bit_list


def paddingToLen(l: list[int], length: int, pos: int = 0):
    while len(l) < length:
        l.insert(pos, 0)
    return l


class HuffmanTable:
    """
    A Huffman Table class
    """

    def __init__(self):
        self.root = []
        self.elements = []

    def BitsFromLengths(self, root, element, pos):
        if isinstance(root, list):
            if pos == 0:
                if len(root) < 2:
                    root.append(element)
                    return True
                return False
            for i in [0, 1]:
                if len(root) == i:
                    root.append([])
                if self.BitsFromLengths(root[i], element, pos - 1) == True:
                    return True
        return False

    def GetHuffmanBits(self, lengths, elements):
        self.elements = elements
        ii = 0
        for i in range(len(lengths)):
            for j in range(lengths[i]):
                self.BitsFromLengths(self.root, elements[ii], i)
                ii += 1

    def Find(self, st):
        r = self.root
        while isinstance(r, list):
            r = r[st.GetBit()]
        return r

    def GetCode(self, st):
        while True:
            res = self.Find(st)
            if res == 0:
                return 0
            elif res != -1:
                return res


class Stream:
    """
    A bit stream class with convenience methods
    """

    def __init__(self, data):
        self.data = data
        self.pos = 0

    def GetBit(self):
        b = self.data[self.pos >> 3]
        s = 7 - (self.pos & 0x7)
        self.pos += 1
        return (b >> s) & 1

    def GetBitN(self, l):
        val = 0
        for _ in range(l):
            val = val * 2 + self.GetBit()
        return val

    def len(self):
        return len(self.data)

    def modifyBit(self, pos: int, bit: int) -> None:
        byte = self.data[pos >> 3]
        bit_list = paddingToLen(intToBitList(byte), 8)
        bit_list[pos & 7] = bit
        self.data[pos >> 3] = bitListToInt(bit_list)


class Decoder:
    def __init__(self, image: str):
        """
        open a JPEG image
        :param image: file name of image to be encrypted
        """
        # open an image file, read data from it
        with open(image, 'rb') as f:
            self.data = list(f.read())

        # initialize variables
        self.huffman_tables: dict[int, HuffmanTable] = {}
        self.height = 0
        self.width = 0
        self.segments = {}

    def separateSegment(self, data: list[int]) -> dict[int, list[int]]:
        """
        separate segments according to markers starting with 0xFF
        """
        # list of markers' indices in self.data
        marker_pos = []
        for i in range(1, len(data)):
            if data[i - 1] == 0xFF and data[i] != 0x00:
                # find a marker
                marker_pos.append(i - 1)
        # list of segments
        segments = {}
        for i in range(1, len(marker_pos)):
            marker = byteListToInt(data[marker_pos[i - 1]: marker_pos[i - 1] + 2])
            segments[marker] = remove0xFF00(data[marker_pos[i - 1]: marker_pos[i]])
        # add the last segment, which starts from the last marker to end
        marker = byteListToInt(data[marker_pos[-1]: marker_pos[-1] + 2])
        segments[marker] = remove0xFF00(data[marker_pos[-1]:])
        return segments

    def decodeHuffman(self, segment: list[int]) -> None:
        """
        segment DHT(Define of Huffman Table) structure:
        -------------------------------------------------------------------
        marker                  2 bytes,
                                0xff, 0xc4 to identify DHT marker
        length                  2 bytes,
                                the length of Huffman table
                                (including length, but not including marker)
        HT information          1 byte,
                                bit 0..3: number of HT (0..3, otherwise error)
                                bit 4: type of HT, 0 = DC table, 1 = AC table
                                bit 5..7: not used, must be 0
        Number of Symbols       16 bytes,
                                Number of symbols with codes of length 1..16,
                                the sum(n) of these bytes is the total number of codes,
                                which must be <= 256
        Symbols                 n bytes,
                                Table containing the symbols in order of increasing code length
                                ( n = total number of codes )
        -------------------------------------------------------------------
        :param segment: DHT segment
        :return:
        """
        marker = byteListToInt(segment[0: 2])
        if marker != 0xFFC4:
            print("Error: invalid DHT segment")
            return
        length = byteListToInt(segment[2: 4])
        data = segment[4: length + 2]
        while len(data) > 0:
            # get information of Huffman table
            offset = 0
            header = data[offset]
            offset += 1

            # get lengths
            lengths = data[offset: offset + 16]
            offset += 16

            # get elements
            elements = []
            for l in lengths:
                elements += data[offset: offset + l]
                offset += l

            hf = HuffmanTable()
            hf.GetHuffmanBits(lengths, elements)
            self.huffman_tables[header] = hf
            data = data[offset:]

    def getImageSize(self, segment: list[int]) -> None:
        """
        Start of Frame structure:
        -------------------------------------------------------------------
        Marker Identifier	    2 bytes, 0xff, 0xc0
        Length	                2 bytes,
        	                    This value equals to 8 + components*3 value
        Data precision	        1 byte,
        	                    This is in bits/sample, usually 8
        	                    (12 and 16 not supported by most software).
        Image height	        2 bytes, this must be > 0
        Image Width         	2 bytes, this must be > 0
        Number of components	1 byte,
        	                    Usually 1 = grey scaled, 3 = color YcbCr or YIQ
        Each component	        3 bytes,
        	                    Read each component data of 3 bytes.
        	                    It contains, (component Id(1byte)(1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q),
        	                    sampling factors (1byte) (bit 0-3 vertical., 4-7 horizontal.),
        	                    quantization table number (1 byte)).
        -------------------------------------------------------------------
        :param segment: SOF segment
        :return:
        """
        marker = byteListToInt(segment[0: 2])
        if marker != 0xFFC0:
            print("Error: invalid SOF segment")
        self.height = byteListToInt(segment[5: 7])
        self.width = byteListToInt(segment[7: 9])

    def decode(self):
        self.segments = self.separateSegment(self.data)
        for marker in self.segments:
            if marker == 0xFFD9:
                # End of Image
                return

            if marker == 0xFFC4:
                # Define of Huffman Table
                self.decodeHuffman(self.segments[marker])
            elif marker == 0xFFC0:
                # Start of Frame0, Baseline DCT
                self.getImageSize(self.segments[marker])


def sortListOnElemLen(list_of_list: list[list[int]]) -> list[int]:
    """
    sort list of list according to their length in descending order
    :param list_of_list: list to be sorted
    :return: sorted index
    """
    len_dict: dict[int, list[int]] = {}
    for i in range(len(list_of_list)):
        elem_len = len(list_of_list[i])
        if elem_len not in len_dict:
            len_dict[elem_len] = []
        len_dict[elem_len].append(i)
    sorted_indices = []
    for i in sorted(len_dict, reverse=True):
        sorted_indices += len_dict[i]
    return sorted_indices


class Crypto:
    def __init__(self, image: str, key: list[int]):
        """
        set value to initialize
        :param image: file name of image to be encrypted
        :param key: secrete key
        """
        # decode basic information
        self.image = image
        self.jpeg_image_decoder = Decoder(image)
        self.jpeg_image_decoder.decode()

        # set key for RC4 encryption & decryption
        self.key = key

    def getAmplitudeBitList(self, st: Stream, idx: int) -> tuple[list[list[int]], list[int]]:
        """
        extract bits list of coefficients from stream
        :param st: bit stream
        :param idx: Huffman table index
        :return: list of bits
        """
        pos_list: list[int] = []
        # DC coefficient
        # size of coded amplitude in bits
        size = self.jpeg_image_decoder.huffman_tables[0 + idx].GetCode(st)
        pos_list.append(st.pos)

        # list of coded amplitude
        amplitude_bits_list: list[list[int]] = [paddingToLen(intToBitList(st.GetBitN(size)), size)]

        l = 1
        while l < 64:
            symbol1 = self.jpeg_image_decoder.huffman_tables[16 + idx].GetCode(st)
            if symbol1 == 0:
                # End of Block
                break

            if symbol1 > 15:
                l += symbol1 >> 4
            size = symbol1 & 0x0F
            pos_list.append(st.pos)
            amplitude_bits_list.append(paddingToLen(intToBitList(st.GetBitN(size)), size))
            if l < 64:
                l += 1
        return amplitude_bits_list, pos_list

    def encryptCoefficientBitList(self, bits_list: list[list[int]]) -> tuple[list[int], list[int]]:
        """
        encrypt bit list denoting coded DCT coefficients
        :param bits_list: coded DCT coefficients
        :return: encrypted bit list
        """
        # sort AC coefficient bit list
        sorted_indices = sortListOnElemLen(bits_list[1:])
        # add in index of DC coefficient
        sorted_indices = [sorted_indices[i] + 1 for i in range(len(sorted_indices))]
        sorted_indices.insert(0, 0)

        # bit list for encryption
        bit_list_for_encryption: list[int] = []
        for i in range(len(sorted_indices)):
            bit_list_for_encryption += bits_list[sorted_indices[i]]

        # encrypt with RC4
        rc4 = RC4(self.key)
        encrypted_byte_list = rc4.encrypt(bitListToByteList(bit_list_for_encryption))
        encrypted_bit_list = byteListToBitList(encrypted_byte_list)
        return encrypted_bit_list, sorted_indices

    def encryptMatrix(self, st: Stream, idx: int):
        # get amplitudes bits
        amplitude_bits_list, amplitude_pos_list = self.getAmplitudeBitList(st, idx)
        # encrypt amplitude bits
        encrypted_bit_list, sorted_indices = self.encryptCoefficientBitList(amplitude_bits_list)
        encrypted_amplitude_bits_list = []

        # divide encrypted bit list according to amplitude length
        for i in range(len(sorted_indices)):
            encrypted_amplitude_bits = []
            for j in range(len(amplitude_bits_list[sorted_indices[i]])):
                encrypted_amplitude_bits.append(encrypted_bit_list.pop(0))
            encrypted_amplitude_bits_list.append(encrypted_amplitude_bits)

        # write back encrypted amplitude into st.data
        for i in range(len(amplitude_pos_list)):
            for j in range(len(encrypted_amplitude_bits_list[i])):
                st.modifyBit(amplitude_pos_list[sorted_indices[i]], encrypted_amplitude_bits_list[i][j])

    def encryptStartOfScan(self, data: list[int]) -> list[int]:
        """
        Start of Scan structure:
        -------------------------------------------------------------------
        Marker Identifier	    2 bytes, 0xff, 0xda
        Length	                2 bytes, 0x0c
        	                    This value equals to 6 + components*3
        Components Number       1 byte, 1-4
        Selector                components*2 bytes, each selector's structure:
                                ------------------------
                                selector    DC, AC table
                                ------------------------
        Spectral select     	1 byte, 0...63
        Successive approx.	    1 byte, 00
        -------------------------------------------------------------------
        :param data: data of Start of Scan segment
        :return:
        """
        marker = byteListToInt(data[0: 2])
        if marker != 0xFFDA:
            print("Error: invalid Start of Scan segment")
            return None

        components = byteListToInt(data[4: 5])
        data_start = 8 + components*2
        st = Stream(data[data_start:])
        for y in range(self.jpeg_image_decoder.height // 8):
            for x in range(self.jpeg_image_decoder.width // 8):
                self.encryptMatrix(st, 0)  # luminance
                self.encryptMatrix(st, 1)  # chrominance r
                self.encryptMatrix(st, 1)  # chrominance b

        # write back encrypted result into data
        for i in range(len(st.data)):
            data[i + data_start] = st.data[i]
        return data

    def encryptImage(self):
        self.encryptStartOfScan(self.jpeg_image_decoder.segments[0xFFDA])

        # save encrypted image
        with open("encrypted_" + self.image, "wb") as f:
            for marker in self.jpeg_image_decoder.segments:
                self.jpeg_image_decoder.segments[marker] = (
                    self.jpeg_image_decoder.segments[marker][0: 2] +
                    add0xFF00(self.jpeg_image_decoder.segments[marker][2:])
                )

                for elem in self.jpeg_image_decoder.segments[marker]:
                    byte = pack("B", elem)
                    f.write(byte)
