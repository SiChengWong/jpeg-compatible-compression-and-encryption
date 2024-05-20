from struct import pack
import rc4
import math


def remove0xFF00(l):
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


def add0xFF00(l):
    l_added_0xFF00 = []
    for i in range(len(l)):
        l_added_0xFF00.append(l[i])
        if l[i] == 0xFF:
            l_added_0xFF00.append(0x00)
    return l_added_0xFF00


def removeMarker(l, marker):
    l_without_marker = []
    marker_len = len(intToBitList(marker)) >> 3
    i = 0
    while i < len(l):
        if i + marker_len < len(l) and byteListToInt(l[i: i + marker_len]) == marker:
            i += marker_len
        else:
            l_without_marker.append(l[i])
            i += 1
    return l_without_marker


def byteListToInt(l):
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


def intToByteList(n: int):
    byte_list = []
    while n > 0:
        byte_list.insert(0, n & 0xFF)
        n = n >> 8
    return byte_list


def intToBitList(n: int):
    bit_list = []
    while n > 0:
        bit_list.insert(0, n & 1)
        n = n >> 1
    return bit_list


def bitListToInt(bit_list):
    n = 0
    for b in bit_list:
        n = (n << 1) + b
    return n


def paddingToLen(l, length, pos=0):
    while len(l) < length:
        l.insert(pos, 0)
    return l


def sortListOnElemLen(list_of_list):
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
    def __init__(self, image):
        """
        open a JPEG image
        :param image: file name of image to be encrypted
        """
        # open an image file, read data from it
        with open(image, 'rb') as f:
            self.data = list(f.read())

        # initialize variables
        self.huffman_tables = {}
        self.huffman_table_mapping = {}
        self.huffman_code_dict = {}

        self.quantization_tables = {}
        self.quantization_mapping = {}

        self.subsample_factor = {}

        self.height = 0
        self.width = 0

        self.segments = {}

    def separateSegment(self, data):
        """
        separate segments according to markers starting with 0xFF
        """
        # list of markers' indices in self.data
        marker_pos = []
        for i in range(1, len(data)):
            if data[i - 1] == 0xFF and data[i] != 0x00:
                # find a marker
                marker_pos.append(i - 1)
        marker_pos.append(len(data))
        # list of segments
        segments = {}
        for i in range(len(marker_pos) - 1):
            marker = byteListToInt(data[marker_pos[i]: marker_pos[i] + 2])
            if marker not in segments:
                segments[marker] = []
            segments[marker].append(remove0xFF00(data[marker_pos[i]: marker_pos[i + 1]]))
        return segments

    def getHuffmanCodeDict(self, root, prev_code, header):
        if isinstance(root, int):
            self.huffman_code_dict[header][root] = prev_code
        else:
            for i in range(len(root)):
                self.getHuffmanCodeDict(root[i], prev_code + [i], header)

    def decodeHuffman(self, segment):
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
        for data in segment:
            data = data[4:]
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

        # generate huffman code table
        for header in self.huffman_tables:
            self.huffman_code_dict[header] = {}
            self.getHuffmanCodeDict(
                self.huffman_tables[header].root,
                [],
                header
            )

    def defineQuantizationTables(self, segment) -> None:
        """
        segment Define Quantization Table structure:
        -------------------------------------------------------------------
        Marker Identifier	    2 bytes, 0xff, 0xdb
        Length                  2 bytes, length of QT
        QT information          1 byte, QT number
        Bytes                   n bytes, gives QT values, n = 64*(precision + 1)
        -------------------------------------------------------------------
        :param segment:
        :return:
        """
        for data in segment:
            data = data[4:]
            while len(data) > 0:
                header = data[0]
                self.quantization_tables[header] = data[1: 1 + 64]
                data = data[65:]

    def baselineDCT(self, segment) -> None:
        """
        Start of Frame structure:
        -------------------------------------------------------------------
        Marker Identifier	    2 bytes, 0xff, 0xc0
        Length	                2 bytes,
        	                    This value equals to 8 + components*3 value
        Data precision	        1 byte,
        	                    bits per pixel, usually 0x08
        Image height	        2 bytes, this must be > 0
        Image Width         	2 bytes, this must be > 0
        Number of components	1 byte,
        	                    Usually 1 = grey scaled, 3 = color YcbCr or YIQ
        Each component	        3 bytes,
        	                    component,
        	                    sampling factor, (bit 0...3, vertical; bit 4...7, horizontal)
        	                    quantization table number
        -------------------------------------------------------------------
        :param segment: SOF segment
        :return:
        """
        for data in segment:
            (
                marker,
                length,
                precision,
                self.height,
                self.width,
                num_of_component
            ) = (
                byteListToInt(data[0: 2]),
                byteListToInt(data[2: 4]),
                data[4],
                byteListToInt(data[5: 7]),
                byteListToInt(data[7: 9]),
                data[9]
            )

            # extract components information
            components = data[10:]
            for i in range(num_of_component):
                # component number
                component = components[i * 3]
                # component subsampling factor
                self.subsample_factor[component] = components[i * 3 + 1]
                # component quantization mapping
                self.quantization_mapping[component] = components[i * 3 + 2]

    def starOfScan(self, segment) -> None:
        """
        Start of Scan structure:
        -------------------------------------------------------------------
        Marker Identifier	    2 bytes, 0xff, 0xda
        Length	                2 bytes, 0x0c
        	                    This value equals to 6 + components*3
        Component Number        1 byte
        Selector                2 bytes for each component,
                                Huffman table for each component's DC, AC coefficient
        Spectral select Start   1 byte, 0x00
        Spectral select End     1 byte, 0x3F
        Successive approx.	    1 byte, 0x00
        -------------------------------------------------------------------
        :param segment: data of Start of Scan segment
        :return:
        """
        for data in segment:
            (
                marker,
                length,
                num_of_component
            ) = (
                byteListToInt(data[0: 2]),
                byteListToInt(data[2: 4]),
                data[4]
            )
            data = data[5:]
            for i in range(num_of_component):
                component = data[i * 2]
                self.huffman_table_mapping[component] = data[i * 2 + 1]


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
                self.baselineDCT(self.segments[marker])
            elif marker == 0xFFDB:
                # Define of Quantization Table
                self.defineQuantizationTables(self.segments[marker])
            elif marker == 0xFFDA:
                # Start of Scan
                self.starOfScan(self.segments[marker])


class Crypto:
    def __init__(self, image, key, iv):
        """
        set value to initialize
        :param iv: initial vector
        :param image: file name of image to be encrypted
        :param key: secrete key
        """
        # decode basic information
        self.image = image
        self.jpeg_image_decoder = Decoder(image)
        self.jpeg_image_decoder.decode()

        self.key = key
        # working mode: CRT
        self.counter = iv

        # compression tag
        self.compression_level = 0

    def getCompressionTag(self):
        marker = 0xFFFE
        if marker in self.jpeg_image_decoder.segments:
            self.compression_level = self.jpeg_image_decoder.segments[marker][0][4]

    def getAmplitudeBitList(self, st, idx):
        """
        extract bits list of coefficients from stream
        :param st: bit stream
        :param idx: Huffman table index
        :return: list of bits
        """
        pos_list = []
        # DC coefficient
        # size of coded amplitude in bits
        size = self.jpeg_image_decoder.huffman_tables[0 + idx].GetCode(st)
        pos_list.append(st.pos)

        # list of coded amplitude
        amplitude_bits_list = [paddingToLen(intToBitList(st.GetBitN(size)), size)]

        l = 1
        while l < 64:
            symbol1 = self.jpeg_image_decoder.huffman_tables[16 + idx].GetCode(st)
            if symbol1 == 0:
                # End of Block
                break

            if symbol1 > 15:
                l += symbol1 >> 4
            size = symbol1 & 0x0F
            if size == 0:
                l += 1
                continue
            pos_list.append(st.pos)
            amplitude_bits_list.append(paddingToLen(intToBitList(st.GetBitN(size)), size))
            if l < 64:
                l += 1
        return amplitude_bits_list, pos_list

    def encryptCoefficientBitList(self, bits_list):
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

        # bits list for encryption
        sorted_bits_list = []
        for i in range(len(sorted_indices)):
            sorted_bits_list.append(bits_list[sorted_indices[i]])

        # encrypt with pseudo-random bit stream
        encrypted_bit_list = []
        # set key for RC4 encryption & decryption
        rand_bit_generator = rc4.RC4RandBitGenerator(
            [(self.key[i] + self.counter) & 0xFF for i in range(len(self.key))]
        )
        self.counter += 1

        encrypted_bit_list = []
        for i in range(len(sorted_bits_list)):
            for j in range(len(sorted_bits_list[i])):
                sorted_bits_list[i][j] = sorted_bits_list[i][j] ^ rand_bit_generator.genRandBit()
            for j in range(self.compression_level):
                rand_bit_generator.genRandBit()
            encrypted_bit_list += sorted_bits_list[i]

        return encrypted_bit_list, sorted_indices

    def encryptMatrix(self, st: Stream, idx: int):
        # get amplitudes bits
        amplitude_bits_list, amplitude_pos_list = self.getAmplitudeBitList(st, idx)
        # encrypt amplitude bits
        encrypted_bit_list, sorted_indices = self.encryptCoefficientBitList(amplitude_bits_list)

        # divide encrypted bit list according to amplitude length
        encrypted_amplitude_bits_list = [[] for i in range(len(sorted_indices))]
        for i in range(len(sorted_indices)):
            for j in range(len(amplitude_bits_list[sorted_indices[i]])):
                encrypted_amplitude_bits_list[sorted_indices[i]].append(encrypted_bit_list.pop(0))

        # write back encrypted amplitude into st.data
        for i in range(len(amplitude_pos_list)):
            for j in range(len(encrypted_amplitude_bits_list[i])):
                st.modifyBit(amplitude_pos_list[i] + j, encrypted_amplitude_bits_list[i][j])

    def encryptStartOfScan(self, data):
        components = byteListToInt(data[4: 5])
        data_start = 8 + components*2
        st = Stream(data[data_start:])

        # num of 8x8 blocks
        sample_sizes = [
            (
                (
                    self.jpeg_image_decoder.subsample_factor[i] >> 4,
                    self.jpeg_image_decoder.subsample_factor[i] & 0x0F,
                    (
                            (self.jpeg_image_decoder.subsample_factor[i] >> 4)
                            * (self.jpeg_image_decoder.subsample_factor[i] & 0x0F)
                    )

                )
            ) for i in self.jpeg_image_decoder.subsample_factor
        ]
        mcu_num = (
                math.ceil(
                    self.jpeg_image_decoder.height / (8 * sample_sizes[0][0])
                ) *
                math.ceil(
                    self.jpeg_image_decoder.width / (8 * sample_sizes[0][1])
                )
        )
        for i in range(mcu_num):
            for j in range(sample_sizes[0][2]):
                self.encryptMatrix(st, self.jpeg_image_decoder.huffman_table_mapping[1] & 0x0F)
            for j in range(sample_sizes[1][2]):
                self.encryptMatrix(st, self.jpeg_image_decoder.huffman_table_mapping[2] & 0x0F)
            for j in range(sample_sizes[2][2]):
                self.encryptMatrix(st, self.jpeg_image_decoder.huffman_table_mapping[3] & 0x0F)

        # write back encrypted result into data
        for i in range(len(st.data)):
            data[i + data_start] = st.data[i]
        return data

    def encryptImage(self, image_output = None):
        self.getCompressionTag()
        for data in self.jpeg_image_decoder.segments[0xFFDA]:
            self.encryptStartOfScan(data)

        # save encrypted image
        if image_output is None:
            image_output = "_" + self.image
        with open(image_output, "wb") as f:
            for marker in self.jpeg_image_decoder.segments:
                for i in range(len(self.jpeg_image_decoder.segments[marker])):
                    self.jpeg_image_decoder.segments[marker][i] = (
                        self.jpeg_image_decoder.segments[marker][i][0: 2] +
                        add0xFF00(self.jpeg_image_decoder.segments[marker][i][2:])
                    )
                for i in range(len(self.jpeg_image_decoder.segments[marker])):
                    for elem in self.jpeg_image_decoder.segments[marker][i]:
                        byte = pack("B", elem)
                        f.write(byte)


class Compression:
    def __init__(self, image: str, compression_level: int = 1):
        """
        :param image: file name of image to be compressed
        :param compression_level: how many bits will be discarded
        """
        self.image = image
        self.jpeg_image_decoder = Decoder(self.image)
        self.jpeg_image_decoder.decode()

        self.compression_level = compression_level

    def addCompressionTag(self):
        """
        comment segment structure:
        -------------------------------------------------------------------
        Marker Identifier	    2 bytes, 0xff, 0xfe
        Length                  2 bytes
        Content                 1 byte, image compressed time(s)
        -------------------------------------------------------------------
        :return:
        """
        marker = 0xFFFE
        length = 0x03
        if marker not in self.jpeg_image_decoder.segments:
            self.jpeg_image_decoder.segments[marker] = [
                    intToByteList(marker) +
                    paddingToLen(intToByteList(length), 2) +
                    [0]
                ]
        # compression time plus compression ratio
        self.jpeg_image_decoder.segments[marker][0][4] += self.compression_level

    def modifyQuantTable(self) -> None:
        marker = 0xFFDB
        # modify quant table segment
        # each element in quantization table doubles self.compression_ratio times
        for k in range(len(self.jpeg_image_decoder.segments[marker])):
            offset = 4
            while offset < len(self.jpeg_image_decoder.segments[marker][k]):
                offset += 1
                for i in range(64):
                    self.jpeg_image_decoder.segments[marker][k][offset + i] = min(
                        (self.jpeg_image_decoder.segments[marker][k][offset + i] << self.compression_level), 0xFF
                    )
                offset += 64
        # modify self.quant
        for header in self.jpeg_image_decoder.quantization_tables:
            for i in range(len(self.jpeg_image_decoder.quantization_tables[header])):
                self.jpeg_image_decoder.quantization_tables[header][i] = min(
                    (self.jpeg_image_decoder.quantization_tables[header][i] << self.compression_level), 0xFF
                )

    def getAmplitudeBitList(self, st, idx):
        # decode a matrix
        amplitude_bits = [[] for i in range(64)]
        # DC coefficient
        size = self.jpeg_image_decoder.huffman_tables[0 + idx].GetCode(st)
        amplitude_bits[0] = paddingToLen(intToBitList(st.GetBitN(size)), size)
        # AC coefficients
        l = 1
        while l < 64:
            symbol_1 = self.jpeg_image_decoder.huffman_tables[16 + idx].GetCode(st)
            if symbol_1 == 0:
                # End of Block
                break

            if symbol_1 > 15:
                l += symbol_1 >> 4
            size = symbol_1 & 0x0F
            amplitude_bits[l] = paddingToLen(intToBitList(st.GetBitN(size)), size)
            if l < 64:
                l += 1
        return amplitude_bits

    def compressCoefficientBits(self, amplitude_bits):
        # compress the matrix
        for i in range(len(amplitude_bits)):
            amplitude_bits[i] = amplitude_bits[i][0: max(len(amplitude_bits[i]) - self.compression_level, 0)]
        return amplitude_bits

    def encodeCoefficientBits(self, idx, amplitude_bits):
        # encode the matrix
        coded_bit_list: list[int] = []
        # DC coefficient
        size = len(amplitude_bits[0])
        coded_bit_list += self.jpeg_image_decoder.huffman_code_dict[0 + idx][size]
        coded_bit_list += amplitude_bits[0]
        # AC coefficient
        l = 1
        run_length = 0
        while l < 64:
            if not amplitude_bits[l]:
                # get run-length
                k = l + 1
                while k < 64:
                    if amplitude_bits[k]:
                        break
                    k += 1

                if k == 64:
                    # End of Block
                    coded_bit_list += self.jpeg_image_decoder.huffman_code_dict[16 + idx][0]
                    break
                else:
                    run_length = k - l
                    while run_length > 0x0F:
                        # 16 continuous zeros
                        coded_bit_list += self.jpeg_image_decoder.huffman_code_dict[16 + idx][0xF0]
                        run_length -= 16
                l = k
            else:
                # non-zero AC coefficient

                # code symbol-1
                symbol_1 = (run_length << 4) + len(amplitude_bits[l])
                coded_bit_list += self.jpeg_image_decoder.huffman_code_dict[16 + idx][symbol_1]
                # code amplitude
                coded_bit_list += amplitude_bits[l]

                # reset run-length
                run_length = 0
                if l < 64:
                    l += 1
        return coded_bit_list

    def compressMatrix(self, st, idx):
        amplitude_bits = self.getAmplitudeBitList(st, idx)
        compressed_amplitude_bits = self.compressCoefficientBits(amplitude_bits)
        coded_compressed_bit_list = self.encodeCoefficientBits(idx, compressed_amplitude_bits)
        return coded_compressed_bit_list

    def compressStartOfScan(self):
        marker = 0xFFDA
        for k in range(len(self.jpeg_image_decoder.segments[marker])):
            data = self.jpeg_image_decoder.segments[marker][k]
            components = byteListToInt(data[4: 5])
            data_start = 8 + components*2
            st = Stream(data[data_start:])
            # num of 8x8 blocks
            sample_sizes = [
                (
                        (
                            self.jpeg_image_decoder.subsample_factor[i] >> 4,
                            self.jpeg_image_decoder.subsample_factor[i] & 0x0F,
                            (
                                    (self.jpeg_image_decoder.subsample_factor[i] >> 4)
                                    * (self.jpeg_image_decoder.subsample_factor[i] & 0x0F)
                            )

                        )
                ) for i in self.jpeg_image_decoder.subsample_factor
            ]
            mcu_num = (
                    math.ceil(
                        self.jpeg_image_decoder.height / (8 * sample_sizes[0][0])
                    ) *
                    math.ceil(
                        self.jpeg_image_decoder.width / (8 * sample_sizes[0][1])
                    )
            )

            # compress each block
            compressed_bit_list = []
            for i in range(mcu_num):
                for j in range(sample_sizes[0][2]):
                    compressed_bit_list += self.compressMatrix(
                        st, self.jpeg_image_decoder.huffman_table_mapping[1] & 0x0F
                    )
                for j in range(sample_sizes[1][2]):
                    compressed_bit_list += self.compressMatrix(
                        st, self.jpeg_image_decoder.huffman_table_mapping[2] & 0x0F
                    )
                for j in range(sample_sizes[2][2]):
                    compressed_bit_list += self.compressMatrix(
                        st, self.jpeg_image_decoder.huffman_table_mapping[3] & 0x0F
                    )

            # write back compressed data
            compressed_bytes = []
            while len(compressed_bit_list) % 8 != 0:
                compressed_bit_list.append(0)
            for i in range(len(compressed_bit_list) // 8):
                compressed_bytes.append(bitListToInt(compressed_bit_list[i * 8: i * 8 + 8]))

            self.jpeg_image_decoder.segments[marker][k] = data[0: data_start] + compressed_bytes

    def compressImage(self, image_output = None):
        self.modifyQuantTable()
        self.compressStartOfScan()
        self.addCompressionTag()

        # save the compressed image
        if image_output is None:
            image_output = "-" + self.image
        with open(image_output, "wb") as f:
            for marker in self.jpeg_image_decoder.segments:
                for i in range(len(self.jpeg_image_decoder.segments[marker])):
                    self.jpeg_image_decoder.segments[marker][i] = (
                        self.jpeg_image_decoder.segments[marker][i][0: 2] +
                        add0xFF00(self.jpeg_image_decoder.segments[marker][i][2:])
                    )
                for i in range(len(self.jpeg_image_decoder.segments[marker])):
                    for elem in self.jpeg_image_decoder.segments[marker][i]:
                        byte = pack("B", elem)
                        f.write(byte)

if __name__ == "__main__":
    filename = 'img/ustc-banner.jpg'
    key = [1,9,5,8,0,9,2,0]
    iv = 2022
    
    encryption = Crypto(filename,key,iv)
    filename_encrypted = 'img/ustc-banner-encrypted.jpg'
    encryption.encryptImage(filename_encrypted)

    decryption = Crypto(filename_encrypted,key,iv)
    filename_decrypted = 'img/ustc-banner-decrypted.jpg'
    decryption.encryptImage(filename_decrypted)
